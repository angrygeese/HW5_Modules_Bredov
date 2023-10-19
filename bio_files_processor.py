import os
import modules.fastq_processor as fp
import modules.bio_files_processor_funcs as bfp_f
from typing import Tuple, Union

def convert_multiline_fasta_to_oneline(input_file: str, output_filename: str = '', output_path: str = 'results', replace_output: bool = True) -> None:
    """Converts multiline fasta sequence to oneline.

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


# умеет доставать n трансляций до гена интереса в грязном виде
def select_genes_from_gbk_to_fasta(input_gbk: str, genes: tuple = (), n_before: int = 3, n_after: int = 1, output_fasta: str = '', output_path: str = ''):
    output_filename = process_paths(input_gbk, output_fasta)
    with open(input_gbk, mode='r') as inp_gbk:
        append_, cache = False, []
        for line in inp_gbk:
            line = line.strip().split()
            if 'CDS' in line:
                if cache:
                    cached_lines = cache[-1][interval]
                    if any (elem for elem in cached_lines if genes in elem):
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

