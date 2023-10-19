import modules.dna_rna_tools as drt
import modules.protein_analyzer_tool as pat
import modules.fastq_processor as fp
from typing import Tuple, Union


def run_dna_rna_tools(*args: str) -> Tuple[list, list]:
    """Provides interface for one from four operations under each given
        nucleotide sequence:
        - get transcribed sequence (`transcribe` operation);
        - get reversed sequence (`reverse` operation);
        - get complementary sequence (`complement` operation);
        - get reversed complementary sequence (`reverse_complement` operation);

    Args:
        *args (str): various number of nucleotide sequences terminated by
            one desired operation.

    Returns tuple containing two list:
        result (list): list with operation results for each valid sequence;
    """
    *seqs, operation = args
    if operation not in drt.OPERATIONS:
        raise ValueError(f'Unknown operation `{operation}`. Please, select from: `transcribe`, `reverse`, `complement`, `reverse_complement`')

    result, corrupt_seqs = [], []
    for seq_index, seq in enumerate(seqs):
        is_seq_valid, complement_dict = drt.check_seq(seq)
        if is_seq_valid:
            result.append(drt.OPERATIONS[operation](seq, complement_dict))
        elif not is_seq_valid:
            corrupt_seqs.append((seq_index, seq))

    drt.print_result(result, corrupt_seqs)
    global TRANSCRIBE_RNA_SEQ
    TRANSCRIBE_RNA_SEQ = False

    result = result[0] if len(result) == 1 else result
    return result


def run_protein_analyzer_tool(*args: str, abbreviation: int = 1) -> Tuple[list, list]:
    """Provides interface for one from five operations under each given aminoacid sequence:
        - get aminoacid content in % (`content_check` operation);
        - get sequence length (`seq_length` operation);
        - get empirical formula (`protein_formula` operation)
        - get protein mass (`protein_mass` operation)
        - get protein charge (`charge` operation)

    Args:
        *args (str): various number of protein sequences terminated by
            one desired operation.
        abbrevition (ште): number of letters in aminocids abbreviation, 1 for
            1-letter and 3 for 3-letter. Defaults to 1.

    Returns:
        tuple(result, corrupt_seqs):
            result (list): list with operation results for each valid sequence;
            corrupt_seqs (list): list with tulpes non-valid sequences and their
                indices;
    """
    *seqs, operation = args
    if operation not in pat.OPERATIONS:
        raise ValueError(f'Unknown operation `{operation}`. Please, select from: `content_check`, `seq_length`, `protein_formula`, `protein_mass`, `charge`')

    result, corrupt_seqs = [], []
    for seq_index, seq in enumerate(seqs):
        is_seq_valid, seq_alt = pat.check_and_procees_seq(seq, abbreviation)
        if is_seq_valid:
            result.append(pat.OPERATIONS[operation](seq_alt))
        elif not is_seq_valid:
            corrupt_seqs.append((seq_index, seq))

    pat.print_result(result, corrupt_seqs)

    result = result[0] if len(result) == 1 else result
    corrupt_seqs = corrupt_seqs[0] if len(corrupt_seqs) == 1 else corrupt_seqs
    return result, corrupt_seqs


def run_fastq_processor(input_path: str, output_filename: str = None, output_path: str = 'fastq_filtrator_results', gc_bounds: Union[Tuple[int, int], int] = (0, 100), length_bounds: Tuple[int, int] = (0, 2**32), quality_thershold: int = 0):
    """Filters reads presented in input fasta file into dictionary using
        three metrics:
        - GC-content;
        - sequence length;
        - average phred quality.
        Then writes filtered reads to file named as specified by
        `output_filename` and stores it in `output_path` folder inside
        directory with input file. If no output filename specified, takes
        name of input file.

    Args:
        input_file (str): path to input file.
        output_filename (str): name of output file. Defaults to None.
        output_path (str): name of folder to store output file. Defaults
            to 'fastq_filtrator_results'.
        gc_bounds (tuple|int): desired interval for GC-content if argument 
            is tuple; upper bound of this interval if argument is integer
            (lower will be 0). Defaults to (0, 100).
        length_bounds (tuple|int): desired interval for sequence length if 
            argument is tuple; upper bound of this interval if argument is 
            integer (lower will be 0). Defaults to (0, 2**32).
        quality_thershold (tuple|int): desired interval for average phred
            quality if argument is tuple; upper bound of this interval if 
            argument is integer (lower will be 0). Defaults to 0.

    Returns:
        None.
    """
    output_filename = fp.process_paths(input_path, output_filename, output_path)
    fastq_dict, result = fp.process_file(input_path), {}
    for name, seq in fastq_dict.items():
        is_seq_valid = fp.check_seq_and_bounds(seq, gc_bounds, length_bounds, quality_thershold)
        if is_seq_valid:
            result[name] = seq
    fp.save_output(result, output_filename)
