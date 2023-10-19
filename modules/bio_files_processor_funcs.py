def write_seqs_dict(output_file: str, seqs_dict: dict, mode: str = 'w') -> None:
    """Consequently writes each dictionary value to file as individual line.
    
    Args:
        output_file (str): path to output file.
        seqs_dict (dict): dictionary containig reads names as keys and tuple
            with reads array and phred quality array.
        mode (str): one-letter-string specifiyng write mode, 'w' for 'write'
            and 'a' for 'append'.
    
    Returns:
        None.
    """
    with open(output_file, mode = mode) as file_write:
        for key in seqs_dict:
            line = seqs_dict[key]
            if line:
                file_write.write(f'{"".join(line)}\n')


def clear_seqs_dict(seqs_dict: dict) -> None:
    """Clears dictionary values if values are lists.
    
    Args:
        output_file (str): path to output file.
        seqs_dict (dict): dictionary containig reads names as keys and tuple
            with reads array and phred quality array.
        mode (str): one-letter-string specifiyng write mode, 'w' for 'write'
            and 'a' for 'append'.
    
    Returns:
        None.
    """
    for key in seqs_dict:
        value = seqs_dict[key]
        if isinstance(value, list):
            value.clear()