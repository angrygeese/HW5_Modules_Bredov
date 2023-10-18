import os
import modules.fastq_processor as fp
import modules.bio_files_processor_funcs as bfp_f
from typing import Tuple, Union

def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None, replace_output=True):
    output_fasta = fp.process_paths(input_fasta, output_fasta, '')
    if replace_output and os.path.isfile(output_fasta):
        os.remove(output_fasta)
    with open(input_fasta, mode='r') as inp_fa:
        seqs_dict = {'name': [], 'seq': []}
        for line in inp_fa:
            line = line.strip()         
            if line.startswith('>'):
                bfp_f.write_seqs_dict(output_fasta, seqs_dict)
                bfp_f.clear_seqs_dict(seqs_dict)
                seqs_dict['name'].append(line)
            # else:
            elif set(line).issubset("ATGCatgc"):
                seqs_dict['seq'].append(line)
        bfp_f.write_seqs_dict(output_fasta, seqs_dict)


def change_fasta_start_pos(input_fasta, shift, output_fasta=None):
    output_fasta = fp.process_paths(input_fasta, output_fasta, '')
    if os.path.isfile(input_fasta):
        fa_dict = {'name': [], 'seq': []}
        with open(input_fasta, mode='r') as inp_fa:
            fa_dict['name'] = inp_fa.readline()
            seq = inp_fa.readline()
            fa_dict['seq'] = seq[shift:] + seq[:shift]
            fp.write_seqs_dict(output_fasta, fa_dict['seq'])
