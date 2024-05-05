from __future__ import annotations


import os
import requests
import sys


from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from collections.abc import Sequence
from datetime import datetime
from dotenv import load_dotenv
from io import StringIO
from typing import Optional, Tuple, Union


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class BiologicalSequenceInitializer(BiologicalSequence):
    def __init__(self):
        self._seq_set = None
        self._seq = None

    # инициализируем алфавит
    @property
    def seq_set(self) -> set:
        return self._seq_set

    @seq_set.setter
    def seq_set(self, seq_set: set):
        if not isinstance(seq_set, set):
            raise TypeError("Must be a set.")
        self._seq_set = seq_set

    # инициализируем последовательность
    @property
    def seq(self) -> str | Sequence[str]:
        return self._seq

    @seq.setter
    def seq(self, seq: str | Sequence[str]):
        if not (
            isinstance(seq, Sequence) and all(isinstance(elem, str) for elem in seq)
        ):
            raise TypeError("Must be a string or sequence of string objects.")
        if not set(seq).issubset(self._seq_set):
            raise TypeError(
                f"Sequence does not match alphabet for {repr(self)}. Please, use only: {self._seq_set}."
            )
        self._seq = seq

    # реализуем заданные в абстрактном классе методы
    def __len__(self):
        return self.seq.__len__()

    def __getitem__(self, index: int | slice):
        return self.seq.__getitem__(index)

    def __str__(self):
        seq = self.seq
        return "".join(seq) if not seq is None else ""

    def __repr__(self):
        return f"BiologicalSequenceInitializer({str(self.seq)})"

    def check_alphabet(self) -> bool:
        return set(self._seq).issubset(self._seq_set)


class NucleicAcidSequence(BiologicalSequenceInitializer):
    def __init__(self):
        super().__init__()  # https://stackoverflow.com/a/64504667
        self._complement_rule = None

    # инициализируем правило комплементарности
    @property
    def complement_rule(self) -> dict:
        return self._complement_rule

    @complement_rule.setter
    def complement_rule(self, complement_rule: dict):
        if not isinstance(complement_rule, dict):
            raise TypeError("Must be a dict.")
        self._complement_rule = complement_rule

    def __repr__(self):
        return f"NucleicAcidSequence({str(self.seq)})"

    def complement(self) -> RNASequence | DNASequence:
        self_class, self_seq_class = type(self), type(self.seq)
        complement_seq = str(self).translate(str.maketrans(self.complement_rule))
        return self_class(self_seq_class(complement_seq))

    def gc_content(self) -> float:
        return sum(map(str(self).count, ("G", "C", "g", "c"))) / len(self)


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq: str | Sequence[str]):
        super().__init__()
        self.seq_set = set("AUGCaugc")
        self.complement_rule = dict(zip("AUGCaugc", "UACGuacg"))
        self.seq = seq

    def __repr__(self):
        return f"RNASequence({str(self.seq)})"


class DNASequence(NucleicAcidSequence):
    def __init__(self, seq: str | Sequence[str]):
        super().__init__()
        self.seq_set = set("ATGCatgc")
        self.complement_rule = dict(zip("ATGCatgc", "TACGtacg"))
        self.seq = seq

    def __repr__(self):
        return f"DNASequence({self.seq})"

    def transcribe(self) -> RNASequence:
        return RNASequence(str(self).translate(str.maketrans(("Tt"), ("Uu"))))


class AminoAcidSequence(BiologicalSequenceInitializer):

    aa_uniprot_content = {
        "A": 9.03,
        "R": 5.84,
        "N": 3.79,
        "D": 5.47,
        "C": 1.29,
        "Q": 3.80,
        "E": 6.24,
        "G": 7.27,
        "H": 2.22,
        "I": 5.53,
        "L": 9.85,
        "K": 4.93,
        "M": 2.33,
        "F": 3.88,
        "P": 4.99,
        "S": 6.82,
        "T": 5.55,
        "W": 1.30,
        "Y": 2.88,
        "V": 6.86,
    }

    def __init__(self, seq: str | Sequence[str]):
        super().__init__()
        self.seq_set = set("ARNDCQEGHILKMFPSTWYV")
        self.seq = seq

    def __repr__(self):
        return f"AminoAcidSequence({str(self.seq)})"

    def content_check(self) -> dict:
        "Returns aminoacids content of the protein"
        seq_content = dict.fromkeys(AminoAcidSequence.aa_uniprot_content.keys(), 0)
        for aacd in str(self.seq).upper():
            seq_content[aacd] = seq_content[aacd] + 1

        seq_length = len(self.seq)
        for aacd, occurence in seq_content.items():
            seq_content[aacd] = 100 * occurence / seq_length

        return seq_content


def gc(sequence):
    return 100 * gc_fraction(sequence, ambiguous='ignore')


def process_paths(
        input_file: str,
        output_filename: str = '',
        output_folder: str = 'results'
        ) -> str:
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


def run_fastq_processor(
        input_file: str,
        *,
        output_filename: Optional[str] = None,
        output_path: str = 'fastq_filtrator_results',
        gc_bounds: Union[Tuple[int, int], int] = (0, 100),
        length_bounds: Tuple[int, int] = (0, 2**32),
        quality_threshold: int = 0
        ):
    """Filters reads presented in input fasta file using three metrics:
    - GC-content;
    - sequence length;
    - average phred quality.

    Args:
        input_file (str): path to input file.
        output_filename (str): name of output file.
        output_path (str): name of folder to store output file.
        gc_bounds (tuple|int): desired interval for GC-content if argument
            is tuple; upper bound of this interval if argument is integer
            (lower will be 0).
        length_bounds (tuple|int): desired interval for sequence length if
            argument is tuple; upper bound of this interval if argument is
            integer (lower will be 0).
        quality_thershold (tuple|int): desired interval for average phred
            quality if argument is tuple; upper bound of this interval if
            argument is integer (lower will be 0).

    Returns:
        None.
    """
    output_filename = process_paths(input_file, output_filename, output_path)
    records, good_reads = SeqIO.parse(input_file, "fastq"), []

    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    min_gc_content, max_gc_content = gc_bounds
    length_min, length_max = length_bounds

    for rec in records:
        rec_len, rec_gc_content = len(rec.seq), gc(rec.seq)
        if (
            min(rec.letter_annotations["phred_quality"]) >= quality_threshold
            and max_gc_content >= rec_gc_content >= min_gc_content
            and length_max >= rec_len >= length_min
        ):
            good_reads.append(rec)
    result_quality = SeqIO.write(good_reads, output_filename, "fastq")
    print(result_quality)


def telegram_logger(chat_id):
    def measure_time(func):
        def wrapper(*args, **kwargs):
            with StringIO() as log_stream:
                redir_stdout, redir_stderr = sys.stdout, sys.stderr
                sys.stderr = sys.stdout = log_stream

                try:
                    env_file_exists, tg_api_token = load_dotenv(), os.getenv('TG_API_TOKEN')

                    if env_file_exists and tg_api_token:
                        start = datetime.now()
                        result = func(*args, **kwargs)
                        end = datetime.now()
                        message = f'Function <code>{func.__name__}</code> sucessfully finished in <code>{end - start}</code>.'
                    else:
                        raise FileNotFoundError(f'''No Telegram API token found in working directory.\
                                                Working directory: <code>{os.getcwd()}</code>\
                                                Enviroment file found: <code>{env_file_exists}</code>\
                                                API token: <code>{tg_api_token}</code>''')
                except Exception as ex_text:
                    message = f'Function <code>{func.__name__}</code> failed with following exeption:\n\n<code>{type(ex_text).__name__}: {ex_text}</code>.'

                sys.stdout, sys.stderr = redir_stdout, redir_stderr

                log_stream_value = log_stream.getvalue()
                if log_stream_value:
                    requests.post(f'https://api.telegram.org/bot{tg_api_token}/sendDocument',
                                params={'chat_id': chat_id, 'caption': message, 'parse_mode': 'HTML'},
                                files={'document': (f'{func.__name__}.log', log_stream_value)}).text
                else:
                    requests.get(f"https://api.telegram.org/bot{tg_api_token}/sendMessage?chat_id={chat_id}&text={message}&parse_mode=HTML").text
                del tg_api_token

                return result
                
        return wrapper
    return measure_time
