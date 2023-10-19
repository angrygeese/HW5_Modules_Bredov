def write_seqs_dict(output_fasta: str, seqs_dict: dict, mode='w'):
    with open(output_fasta, mode=mode) as file_write:
        for name in seqs_dict:
            line = seqs_dict[name]
            if line:
                file_write.write(f'{"".join(line)}\n')


def clear_seqs_dict(seqs_dict: dict):
    seqs_dict['name'].clear()
    seqs_dict['seq'].clear()