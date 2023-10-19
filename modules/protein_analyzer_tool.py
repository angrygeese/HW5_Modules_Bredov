from typing import Iterable, Tuple

AA_TR_DICT = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
    'Cys': 'C', 'Gln': 'E', 'Glu': 'Q', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
AA_UNIPROT_CONTENT = {
    'A': 9.03, 'R': 5.84, 'N': 3.79, 'D': 5.47, 'C': 1.29,
    'Q': 3.80, 'E': 6.24, 'G': 7.27, 'H': 2.22, 'I': 5.53,
    'L': 9.85, 'K': 4.93, 'M': 2.33, 'F': 3.88, 'P': 4.99,
    'S': 6.82, 'T': 5.55, 'W': 1.30, 'Y': 2.88, 'V': 6.86
    }
AA_MOL_MASS = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
    'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
    'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
    'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
AA_FORMULA = {
    'A': (3, 7, 1, 2, 0), 'R': (6, 14, 4, 2, 0), 'N': (4, 8, 2, 3, 0), 'D': (4, 7, 1, 4, 0), 'C': (3, 7, 1, 2, 1),
    'Q': (5, 10, 2, 3, 0), 'E': (5, 9, 1, 4, 0), 'G': (2, 5, 1, 2, 0), 'H': (6, 9, 3, 2, 0), 'I': (6, 13, 1, 2, 0),
    'L': (6, 13, 1, 2, 0), 'K': (6, 14, 2, 2, 0), 'M': (5, 11, 1, 2, 1), 'F': (9, 11, 1, 2, 0), 'P': (5, 9, 1, 2, 0),
    'S': (3, 7, 1, 3, 0), 'T': (4, 9, 1, 3, 0), 'W': (11, 12, 2, 2, 0), 'Y': (9, 11, 1, 3, 0), 'V': (5, 11, 1, 2, 0)
    }
AA_CHARGES = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
    'Q': 0, 'E': -1, 'G': 0, 'H': 1, 'I': 0,
    'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
    }


def mann_whitney_u(seq1: Iterable[int], seq2: Iterable[int]) -> bool:
    """
    Mann-Whitney U-test. Used to compare aminoacids composition in given sequence with average protein composition provided by Uniprot.
    Used as a second step of checkup in `check_and_procees_seq` function if sequence uses 1-letter abbreviation.
    """
    len_seq1, len_seq2 = len(seq1), len(seq2)
    r1, r2 = dict.fromkeys(map(str, seq1), 0), dict.fromkeys(map(str, seq2), 0)

    r = sorted(list(seq1) + list(seq2))
    r_dict = dict.fromkeys(map(str, r), ())

    for index, value in enumerate(r):
        value = str(value)
        r_dict[value] = r_dict[value] + (index + 1,)
    for elem in r_dict:
        r_dict[elem] = sum(r_dict[elem]) / len(r_dict[elem])

    for value in seq1:
        value = str(value)
        r1[value] = r1[value] + r_dict[value]

    for value in seq2:
        value = str(value)
        r2[value] = r2[value] + r_dict[value]

    u1 = (len_seq1 * len_seq2) + len_seq1 * (len_seq1 + 1) / 2 - sum(r1.values())
    u2 = (len_seq1 * len_seq2) + len_seq2 * (len_seq2 + 1) / 2 - sum(r2.values())

    u_stat = min(u1, u2)

    if u_stat <= 127:
        return False

    return True


def decomposition(seq: str) -> list:
    """
    Decomposes 3 letter-sequence into list of single aminoacids.
    Used in `check_and_procees_seq` function if sequence uses 3-letter abbreviation.
    """
    len_seq, dec_seq = len(seq), []
    if len(seq) % 3 == 0:
        for i in range(0, len_seq, 3):
            dec_seq.append(seq[i:i+3].lower().capitalize())

    return dec_seq


def seq_transform(seq: list) -> str:
    """
    Transforms aminoacid abbreviation format from 3-letter into 1-letter.
    Used in `check_and_procees_seq` function if sequence uses 3-letter abbreviation.
    """
    seq_tr = ''
    for aa in seq:
        seq_tr += AA_TR_DICT[aa]

    return seq_tr


def aa_content_check(seq: str) -> dict:
    "Returns aminoacids content of the protein"
    seq_content = dict.fromkeys(AA_UNIPROT_CONTENT.keys(), 0)
    for aacd in seq.upper():
        seq_content[aacd] = seq_content[aacd] + 1

    seq_length = len(seq)
    for aacd, occurence in seq_content.items():
        seq_content[aacd] = 100 * occurence / seq_length

    return seq_content


def check_and_procees_seq(seq: str, abbreviation: int = 1) -> Tuple[bool, str]:
    "Checks whether the string is valid protein sequence and transforms to 1-letter abbreviation if required"
    exit_code = False
    if isinstance(seq, str):
        if abbreviation == 3:
            seq = decomposition(seq)
            exit_code = bool(seq) and set(seq).issubset(set(AA_TR_DICT.keys()))
            if exit_code:
                seq = seq_transform(seq)
        elif abbreviation == 1:
            seq_set = set(seq.upper())
            exit_code = bool(seq) and seq_set.issubset(set(AA_TR_DICT.values()))
            if exit_code:
                seq_content, uniprot_content = aa_content_check(seq).values(), AA_UNIPROT_CONTENT.values()
                is_seq_content_valid = mann_whitney_u(seq_content, uniprot_content) if len(seq_set) == 20 else True
                exit_code = is_seq_content_valid
        else:
            raise ValueError('Incorrect abbreviation. Must be 1 or 3, type `int`')

    return exit_code, seq


def seq_length(seq: str) -> int:
    "Counts length of the protein"
    return len(seq)


def protein_mass(seq: str) -> float:
    "Counts molecular weight of the protein"
    seq_weight = 0
    for aa in seq:
        seq_weight += AA_MOL_MASS[aa]
    return float(f'{seq_weight:.02f}')


def protein_formula(seq: str) -> dict:
    "Returns molecular formula of the protein"
    seq_formula = (0, 0, 0, 0, 0)
    for aa in seq:
        seq_formula = tuple(map(sum, zip(seq_formula, AA_FORMULA[aa])))
    aa_formula = {atom: number for atom, number in zip('CHNOS', seq_formula)}
    return aa_formula


def aa_chain_charge(seq: str, aa_charges: dict = AA_CHARGES) -> int:
    "Returns charge of the protein (pH=7)"
    aa_charge = 0
    for aa in seq.upper():
        aa_charge += aa_charges[aa]

    return aa_charge


def print_result(result: list, corrupt_seqs: list):
    "Provides brief description for `run_protein_analyzer_tool` output."
    len_seq, len_corr_seq = len(result), len(corrupt_seqs)
    len_seqs = len_seq + len_corr_seq
    success = ["+" for _ in range(len_seqs)]
    if not len_corr_seq:
        print(f'All {len_seqs} sequence(s) processed successfully')
    elif len_corr_seq:
        for i in corrupt_seqs:
            success[i[0]] = "-"
        print(f'Processing result: [{"".join(success)}]\n')
        print(f'{len_seq} sequence(s) out of {len_seq + len_corr_seq} given have been processed successfully.')
        print(f'{len_corr_seq} has been recognized as corrupted, i.e. non-protein')


OPERATIONS = {'content_check': aa_content_check, 'seq_length': seq_length, 'protein_formula': protein_formula, 'protein_mass': protein_mass, 'charge': aa_chain_charge}
