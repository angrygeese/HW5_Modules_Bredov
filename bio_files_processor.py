import os


from _io import TextIOWrapper
from dataclasses import dataclass
from typing import Union


def process_paths(input_file: str, output_filename: str='', output_folder: str='results') -> str:
    """Checks input path and generates name and folder for  input file.

    Args:
        input_file (str): path to input file.
        output_filename (str): name of output file.
        output_folder (str): name of output folder.

    Returns:
        output_file (str): path to output file.
    """
    if os.path.isfile(input_file):
        inp_path, inp_filename = os.path.split(input_file)
        if output_filename and not output_filename.endswith('.fasta'):
            output_filename += '.fasta'
        else:
            output_filename = inp_filename
        output_path = os.path.join(inp_path, output_folder)

        if not os.path.exists(output_path):
            os.mkdir(output_path) 
        output_file = os.path.join(output_path, output_filename)
    else:
        raise FileNotFoundError("Incorrect input path: no such file or directory")

    return output_file


def write_seqs_dict(output_file: str, seqs_dict: dict, mode: str = 'w') -> None:
    """Consequently writes each dictionary value to file as individual line.
    
    Args:
        output_file (str): path to output file.
        seqs_dict (dict): dictionary containig reads names as keys and tuple
            with reads array and phred quality array.
        mode (str): one-letter-string specifiyng write mode, 'w' for 'write'
            and 'a' for 'append'.
    
    Returns:
        None.
    """
    with open(output_file, mode='w') as file_write:
        for name, seq in seqs_dict.items():
            file_write.write(f'{name}\n')
            file_write.write(f'{"".join(seq)}\n')


def convert_multiline_fasta_to_oneline(
        input_file: str,
        output_filename: str = '',
        output_path: str = 'results'
        ) -> None:
    """Converts multiline fasta sequence to oneline.

    Args:
        input_file (str): path to input file.
        output_filename (str): name of output file. Defaults to empty string.
        output_path (str): name of folder to store output file. Defaults to
            'results'. 

    Returns:
        None.
    """     
    output_filename = process_paths(input_file, output_filename, output_path)
    with open(input_file, mode='r') as inp_fa:
        seqs_dict = {}
        for line in inp_fa:
            line = line.strip()         
            if line.startswith('>'):
                seq_id = line
                seqs_dict[seq_id] = []
            else:
                seqs_dict[seq_id].append(line)
            
    write_seqs_dict(output_filename, seqs_dict)


def change_fasta_start_pos(
        input_file: str,
        shift: int = 0,
        output_filename: str = '',
        output_path: str = 'results'
        ) -> None:
    """Shifts starting nucleotide position to desired position for one line fasta file.

    Args:
        input_file (str): path to input file.
        shift (int): new position for starting nucleotide. Defaults to zero.
        output_filename (str): name of output file. Defaults to empty string.
        output_path (str): name of folder to store output file. Defaults to
            'results'. 

    Returns:
        None.
    """   
    output_filename = process_paths(input_file, output_filename, output_path)
    if os.path.isfile(input_file):
        fa_dict = {}
        with open(input_file, mode='r') as inp_fa:
            name = inp_fa.readline().strip()
            seq = inp_fa.readline().strip()
            fa_dict[name] = seq[shift:] + seq[:shift]
            write_seqs_dict(output_filename, fa_dict)


def select_genes_from_gbk_to_fasta(
        input_gbk: str,
        genes: tuple = (),
        n_before: int = 3,
        n_after: int = 1,
        output_fasta: str = '',
        output_path: str = ''):
    """Умеет доставать n трансляций до гена интереса в грязном виде
    """
    output_filename = process_paths(input_gbk, output_fasta)
    with open(input_gbk, mode='r') as inp_gbk:
        append_, cache = False, []
        for line in inp_gbk:
            line = line.strip().split()
            if 'CDS' in line:
                if cache:
                    cached_lines = cache[-1][interval]
                    if any(elem for elem in cached_lines if genes in elem):
                        for CDS in cache:
                            for name in CDS:
                                cached_lines = CDS[name]
                                trans_start = cached_lines.index(*[elem for elem in cached_lines if 'translation' in elem])
                                translation = cached_lines[trans_start:]
                                cached_lines.clear()
                                cached_lines.append(''.join(translation))
                        break
                interval = line[1]
                if len(cache) == n_before:
                    del cache[0]
                cache.append({interval: []})
                append_ = True
                continue
            if 'ORIGIN' in line:
                append_ = False
                if cache:
                    cached_lines = cache[-1][interval]
                    gene = [elem for elem in cached_lines if genes in elem]
                    if gene:
                        break
            if append_:
                cache[-1][interval].extend(line)
        return cache


@dataclass
class FastaRecord:
    rec_id: str
    rec_desc: str
    rec_seq: str

    def __bool__(self):
        if self.rec_id or self.rec_desc:
            return True
        return False

    def __repr__(self):
        return f"ID: {self.rec_id}\nDescription: {self.rec_desc}\nSequence: {self.rec_seq}\n"


class OpenFasta:

    nuc_alphabet = set("ATUGCNatugcn")

    def __init__(self, file: Union[str, os.PathLike], mode: str = "r"):
        self.file = file
        self.mode = mode
        self.handler = None
        self.current = None

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        if isinstance(file, (str, os.PathLike)):
            self._file = file
        else:
            raise ValueError(f"Expected str or os.PathLike object, not {type(file)}")

    def __enter__(self):
        self.file_open()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.file_close()

    def __iter__(self):
        self.file_open()
        return self

    def __next__(self):
        record: FastaRecord = FastaRecord("", "", "")  # `record` хранит `FastaRecord`, полученный в текущей итерации `__next__`
        # логика такая: если текущая строка существует, входим в `while` и итерируемся по `self.handler`
        if self.current is None:
            raise TypeError("I/O operation not possible, open file first.")
        while isinstance(self.current, str):
            # если в строке есть ID или описание (строка начинается с '>'), то проверяем:
            if self.current.startswith(">"):
                # 1) мы дошли до неё в этом вызове `__next__` (`record` не пустой) - завершаем формирование `FastaRecord` и возвращаем его.
                if record:
                    record.rec_seq = "".join(record.rec_seq)
                    return record
                # 2) мы дошли до неё в предыдущем вызове `__next__` (`record` пустой) - создаём новый `FastaRecord`
                rec_id, rec_desc = self.current.split(" ", 1)
                record = FastaRecord(rec_id.lstrip(">"), rec_desc, [])
            # если в строке есть нуклеотиды - пишем их в `rec_seq` текущей `FastaRecord`
            elif set(self.current).issubset(OpenFasta.nuc_alphabet):
                record.rec_seq.append(self.current)
            try:
                self.current = next(self.handler).strip()
            # если в ходе итерации возникает `StopIteration` - перезаписываем текущую строку как `-1`, завершаем формирование `FastaRecord` и возвращаем его.
            except StopIteration:
                self.current = -1
                record.rec_seq = "".join(record.rec_seq)
                return record
        # если текущая строка `-1` - то в цикле выше мы дошли до конца файла, поэтому при дальнейших попытках вызвать `__next__` вызываем `StopIteration`
        if self.current == -1:
            raise StopIteration

    def file_open(self):
        if not isinstance(self.handler, TextIOWrapper) or self.handler.closed:
            self.handler = open(self.file, mode=self.mode)
            self.current = next(self.handler).strip()

    def file_close(self):
        if isinstance(self.handler, TextIOWrapper) and not self.handler.closed:
            self.handler.close()

    def read_record(self):
        return next(self)

    def read_records(self):
        records = []
        next_record_exists = True
        while next_record_exists:
            try:
                next_record = next(self)
            except StopIteration:
                next_record_exists = False
            else:
                records.append(next_record)
        return records
