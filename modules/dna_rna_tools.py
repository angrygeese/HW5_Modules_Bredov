# PR comment
from typing import Tuple, Union

DNA_SET, RNA_SET, NONSENCE_SET = frozenset('ATGC'), frozenset('AUGC'), frozenset('UT')
COMPLEMENT_RNA = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                  'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
COMPLEMENT_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
TRANSCRIBE_RNA_SEQ = False


def check_seq(seq: str, check_set: frozenset = NONSENCE_SET) -> Union[Tuple[bool, dict], None]:
    """Checks whether the string is valid aminoacid sequence.
    
    Args:
        seq (str): sequence to be checked.

    Returns:
        exit_code (bool): returns `True` if sequence is valid, otherwise 
            returns `False`.
        complement_dict (dict|None): dictionary cointaining complementary rule
            for sequence. `None` for non-valid sequences, `COMPLEMENT_DNA` for
            DNA sequence and `COMPLEMENT_RNA` for RNA sequence
    """
    exit_code, complement_dict = False, None
    if isinstance(seq, str):
        exit_code = bool(seq) and not check_set.issubset(seq)
        if exit_code:
            if set(seq.upper()).issubset(DNA_SET):
                complement_dict = COMPLEMENT_DNA
            elif set(seq.upper()).issubset(RNA_SET):
                complement_dict = COMPLEMENT_RNA
    return exit_code, complement_dict


def transcribe(seq: str, complement_dict: dict) -> str:
    """Transcribes sequence.
    
    Args:
        seq (str): sequence to be transcribed.
        complement_dict (dict): dictionary cointaining complementary rule
            for sequence. Used to determine whether sequence is DNA or RNA. 

    Returns:
        transc_seq (str): transcribed sequence; if sequence is RNA, returns
            sequence itself
    """
    global TRANSCRIBE_RNA_SEQ
    if complement_dict is COMPLEMENT_DNA:
        transc_seq = seq.replace('T', 'U').replace('t', 'u')
    elif complement_dict is COMPLEMENT_RNA:
        transc_seq = seq
        if not TRANSCRIBE_RNA_SEQ:
            print('Note: you have passed RNA sequence(s) for transcription. Output will contain RNA sequences themselves.\n')
            TRANSCRIBE_RNA_SEQ = True
    return transc_seq


def reverse(seq: str, complement_dict: dict = None) -> str:
    """Makes reverse sequence.
    
    Args:
        seq (str): sequence to reverse.
        complement_dict (dict): dictionary cointaining complementary rule
            for sequence. Added to provide call unification inside primary
            function, `run_dna_rna_tools` (i.e. all functions were able
            to recieve 2 arguments).

    Returns:
        seq_tr (str): complement sequence.
    """
    return seq[::-1]


def complement(seq: str, complement_dict: dict) -> str:
    """Reconstructs complementary sequence using complementary rule.
    
    Args:
        seq (str): sequence to make complement.
        complement_dict (dict): dictionary cointaining complementary rule
            for sequence.

    Returns:
        seq_tr (str): complementary sequence.
    """
    seq_tr = ''.join([complement_dict[nuc] for nuc in seq])
    return seq_tr


def reverse_complement(seq: str, complement_dict: dict) -> str:
    """Reconstructs reverse complement sequence using complementary rule.
    
    Args:
        seq (str): sequence to make complement.
        complement_dict (dict): dictionary cointaining complementary rule
            for sequence.

    Returns:
        seq_tr (str): reverse complement sequence.
    """
    return reverse(complement(seq, complement_dict))


def print_result(result: list, corrupt_seqs: list):
    "Provides brief description for `run_dna_rna_tools` output."
    len_seq, len_corr_seq = len(result), len(corrupt_seqs)
    len_seqs = len_seq + len_corr_seq
    success = ['+' for _ in range(len_seqs)]
    if not len_corr_seq:
        print(f'All {len_seqs} sequence(s) processed successfully')
    elif len_corr_seq:
        for i in corrupt_seqs:
            success[i[0]] = "-"
        print(f'Processing result: [{"".join(success)}]\n')
        print(f'{len_seq} sequence(s) out of {len_seq + len_corr_seq} given have been processed successfully.')
        print(f'{len_corr_seq} has been recognized as corrupted, i.e. non-DNA as well as non-RNA.')


OPERATIONS = {'transcribe': transcribe, 'reverse': reverse, 'complement': complement, 'reverse_complement': reverse_complement}
