def read_fasta(filepath : str) -> dict:
    """
    This function reads a FASTA file and maps it in a dictionary with protein IDs and their sequence.

    Parameters
    ----------
    filepath: str
        The path to the file
    
    Returns
    ----------
    dict
        Dictionary with Protein IDs as keays and sequences as values (entries).
    """
    protein_map = {}                            # defines empty dictionary 
    current_id = None                           # defines id with no value in beginning
    current_sequence = []                       # defines empty list

    with open(filepath, 'r') as fasta_handle:                             # open file as fasta_handle
        for line in fasta_handle:                                           # goes through all lines
            stripped = line.strip()                                         # removes any leading, and trailing whitespaces
            if not stripped:                                                # if its empty
                continue                                                    # line is skipped
            if stripped.startswith('>'):                                    # if line is a header
                if current_id is not None:                                  # and there was a header before
                    protein_map[current_id] = ''.join(current_sequence)     # adds the sequence to the current dictionary entry in protein_map
                    current_sequence = []                                   # current sequence emptied (now new protein)
                current_id = stripped.split('|')[1]                         # second part of header (between '|') taken as Protein ID
            else:                                                           # if sequence line:
                current_sequence.append(stripped)                           # sequence is appended
 
    if current_id is not None:                                          # if there was an entry before this one:
        protein_map[current_id] = ''.join(current_sequence)             # adds the sequence to the current dictionary entry in protein_map
        
    return protein_map