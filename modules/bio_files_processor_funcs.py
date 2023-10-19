def write_seqs_dict(output_filename: str, seqs_dict: dict, mode: str = 'w') -> None:
    "Consequently writes each dictionary value to file as individual line"
    with open(output_filename, mode = mode) as file_write:
        for key in seqs_dict:
            line = seqs_dict[key]
            if line:
                file_write.write(f'{"".join(line)}\n')


def clear_seqs_dict(seqs_dict: dict) -> None:
    "Clears dictionary values if values are lists"
    for key in seqs_dict:
        value = seqs_dict[key]
        if isinstance(value, list):
            value.clear()