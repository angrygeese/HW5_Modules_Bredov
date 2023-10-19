import os
from typing import Tuple, Union
if __name__ == '__main__':
    from dna_rna_tools import check_seq, COMPLEMENT_DNA
else:
    from .dna_rna_tools import check_seq, COMPLEMENT_DNA


def process_paths(input_file: str, output_filename: str = '', output_folder: str = 'results') -> str:
    """Checks input path and generates name and folder for output file.

    Args:
        input_file (str): path to input file.
        output_filename (str): name of output file. Defaults to empty string.
        output_folder (str): name of output folder. Defaults to 'results'.

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


def process_file(input_file: str, fastq_dict: dict = None) -> dict:
    """Reads fasta file and transforms its content into dictionary.

    Args:
        input_file (str): path to input file.
        fastq_dict (dict): dictionary to store content from input file. Defaults
            to None.

    Returns:
        fastq_dict (dict): dictionary containig reads names as keys and tuple
            with reads array and phred quality array.
    """
    if fastq_dict is None:
        fastq_dict = {}
    with open(input_file, mode='r') as file_read:
        cache = []
        for index, line in enumerate(file_read):
            cache.append(line.strip())
            if (index + 1) % 4 == 0:
                name, seq, comm, phred = cache
                fastq_dict[name] = (seq, comm, phred)
                cache.clear()

    return fastq_dict


def save_output(fastq_dict: dict, output_file: str) -> None:
    """Writes dictionary content into file.

    Args:
        fastq_dict (dict): dictionary containig reads names as keys and tuple
            with reads array and phred quality array as values.
        output_file (str): path to output file.

    Returns:
        None.
    """
    with open(output_file, mode='w') as file_write:
        for name, (nuc_seq, comm, phred_seq) in fastq_dict.items():
            file_write.write(f'{name}\n{nuc_seq}\n{comm}\n{phred_seq}\n')


def is_in_range(value, val_range: Union[Tuple[int, int], int]) -> bool:
    """Checks whether value belongs to desired interval.

    Args:
        value (int): value to be checked.
        val_range (tuple|int): desired interval if argument is tuple;
            upper bound of interval if argument is integer (lower will be 0).

    Returns:
        is_val_in_range (bool): returns `True` if value belongs to desired
            interval, otherwise returns `False`.
    """
    if isinstance(val_range, int):
        val_range = (0, val_range)
    lower_bound, upper_bound = val_range
    is_val_in_range = lower_bound <= value <= upper_bound
    return is_val_in_range


def check_seq_and_bounds(seqs: Tuple[str, str, str], gc_bounds, length_bounds, quality_thershold) -> bool:
    """Counts GC-content and average quality, then checks whether both of them
        belong to specified intervals. 

    Args:
        seqs (tuple): tuple containing three strings, DNA sequence, line
            with comment and phred quality score for each position.
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
        exit_code (bool): returns `True` if:
            - all three metrics - sequence length, GC-content and average
            quality belong to specified interval, othwerwise returns `False`.
            - sequence is empty.
    """
    nuc_seq, _, phred_seq = seqs
    exit_code, complement_dict = check_seq(nuc_seq)
    if exit_code and complement_dict == COMPLEMENT_DNA:
        gc_count, phred_count = 0, 0
        for index, (nuc, phred) in enumerate(zip(nuc_seq, phred_seq)):
            gc_count += nuc in {"G", "C"}
            phred_count += ord(phred) - 33
        if not is_in_range(index + 1, length_bounds):
            return False
        gc_content, quality = round(gc_count * 100 / (index + 1), 2), round(phred_count / (index + 1), 2)
        exit_code = is_in_range(gc_content, gc_bounds) and not is_in_range(quality, quality_thershold)
    return exit_code if nuc_seq else True
