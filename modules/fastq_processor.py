# PR comment
from typing import Tuple, Union
from dna_rna_tools import check_seq, COMPLEMENT_DNA


def is_in_range(value, val_range: Union[Tuple[int, int], int]):
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


def check_seq_and_bounds(seq_pair: Tuple[str, str], gc_bounds, length_bounds, quality_thershold):
    """Counts GC-content and average quality, 
        Then checks whether both of them belong to specified intervals. 

    Args:
        seq_pair (tuple): tuple containing two strings, DNA sequence and 
            hred quality score for each position.
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
        exit_code (bool): returns `True` if all three metrics - sequence
            length, GC-content and average quality belong to specified 
            interval, othwerwise returns `False`.
    """
    nuc_seq, phred_seq = seq_pair
    exit_code, complement_dict = check_seq(nuc_seq)
    if exit_code and complement_dict == COMPLEMENT_DNA:
        gc_count, phred_count = 0, 0
        for index, nuc_and_phred in enumerate(zip(nuc_seq, phred_seq)):
            nuc, phred = nuc_and_phred
            gc_count += nuc in {"G", "C"}
            phred_count += ord(phred) - 33
            if not is_in_range(index + 1, length_bounds):
                return False
        gc_content, quality = round(gc_count * 100 / (index + 1), 2), round(phred_count / (index + 1), 2)
        exit_code = is_in_range(gc_content, gc_bounds) and not is_in_range(quality, quality_thershold)
    return exit_code if nuc_seq else True
