import os
import modules.fastq_processor as fp
import modules.bio_files_processor_funcs as bfp_f
from typing import Tuple, Union

def convert_multiline_fasta_to_oneline(input_file: str, output_filename: str = '', output_path: str = 'results', replace_output: bool = True) -> None:
    """Converts multiline fasta to oneline.

    Args:
        input_file (str): path to input file.
        output_filename (str): name of output file. Defaults to empty string.
        output_path (str): name of folder to store output file. Defaults to
            'results'. 
        replace_output (bool): True to replace output file if it exists, False
            to append output to existing file.

    Returns:
        None.
    """     
    output_filename = fp.process_paths(input_file, output_filename, output_path)
    if replace_output and os.path.isfile(output_filename):
        bfp_f.write_seqs_dict(output_filename, {})
    with open(input_file, mode='r') as inp_fa:
        seqs_dict = {'name': [], 'seq': []}
        for line in inp_fa:
            line = line.strip()         
            if line.startswith('>'):
                bfp_f.write_seqs_dict(output_filename, seqs_dict, 'a')
                bfp_f.clear_seqs_dict(seqs_dict)
                seqs_dict['name'].append(line)
            else:
                seqs_dict['seq'].append(line)
            
        bfp_f.write_seqs_dict(output_filename, seqs_dict, 'a')


def change_fasta_start_pos(input_file: str, shift: int = 0, output_filename: str = '', output_path: str = 'results') -> None:
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
    output_filename = fp.process_paths(input_file, output_filename, output_path)
    if os.path.isfile(input_file):
        fa_dict = {'name': [], 'seq': []}
        with open(input_file, mode='r') as inp_fa:
            fa_dict['name'] = inp_fa.readline().strip()
            seq = inp_fa.readline().strip()
            fa_dict['seq'] = seq[shift:] + seq[:shift]
            bfp_f.write_seqs_dict(output_filename, fa_dict)
