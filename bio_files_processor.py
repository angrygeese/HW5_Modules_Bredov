import os
import modules.fastq_processor as fp
import modules.bio_files_processor_funcs as bfp_f
from typing import Tuple, Union

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None, replace_output: bool = True):
    output_fasta = fp.process_paths(input_fasta, output_fasta, '')
    if replace_output and os.path.isfile(output_fasta):
        bfp_f.write_seqs_dict(output_fasta, {})
    with open(input_fasta, mode='r') as inp_fa:
        seqs_dict = {'name': [], 'seq': []}
        for line in inp_fa:
            line = line.strip()         
            if line.startswith('>'):
                bfp_f.write_seqs_dict(output_fasta, seqs_dict, 'a')
                bfp_f.clear_seqs_dict(seqs_dict)
                seqs_dict['name'].append(line)
            elif set(line).issubset("ATGCatgc"):
                seqs_dict['seq'].append(line)
        bfp_f.write_seqs_dict(output_fasta, seqs_dict, 'a')


def change_fasta_start_pos(input_fasta: str, shift: str, output_fasta: str = None):
    output_fasta = fp.process_paths(input_fasta, output_fasta, '')
    if os.path.isfile(input_fasta):
        fa_dict = {'name': [], 'seq': []}
        with open(input_fasta, mode='r') as inp_fa:
            fa_dict['name'] = inp_fa.readline().strip()
            seq = inp_fa.readline().strip()
            fa_dict['seq'] = seq[shift:] + seq[:shift]
            bfp_f.write_seqs_dict(output_fasta, fa_dict)
