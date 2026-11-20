import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    This function cleaves protein sequences according to cleavage patterns of enzymes.
    Additionally, it filters for peptides in a specific length range.
    
    Parameters
    ----------
    protein_seq: str
        sequence of the protein
    cleave_pattern: re
        pattern after/before which aminocids the enzymes cut
    min_pep_len: int
        minimal length of peptides that will be selected
    max_pep_len: int
            maximal length of peptides that will be selected
    Returns
    -------
    list, peptide sequences
    int, number of peptides
    """
    
    peptides = re.split(cleave_pattern, protein_seq)
    peptides_filt = [len_pep for len_pep in peptides if len(len_pep) in range (min_pep_len, max_pep_len+1) ] # list comprehension for filtering, range includes first but not last value in the brackets
    length = len(peptides_filt)
    #print(peptides_filt)
    #print(f'Nr. of digested peptides: {length}')
    
    return(peptides_filt, length)

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Add a short description here.

    Parameters
    ----------
    protein_map: dict, Dictionary of Protein IDs and according sequences.
    cleave_pattern: re, pattern after/before which aminocids the enzymes cut
    min_pep_len: int, minimal length of peptides that will be selected
    max_pep_len: int, maximal length of peptides that will be selected

    Returns
    -------  
    dict, protein_id : list of peptides obtaied from sequence when cutting
    """
    peptide_dict = {}
    
    for protein_id, protein_seq in protein_map.items():
        peptides = re.split(cleave_pattern, protein_seq)
        peptides_filt = [len_pep for len_pep in peptides if len(len_pep) in range (min_pep_len, max_pep_len+1) ] # list comprehension for filtering, range includes first but not last value in the brackets
        peptide_dict[protein_id] = peptides_filt
           
    return(peptide_dict)


def compute_sequence_coverage(protein_seq, peptides):
    """
    Calculates sequence coverage = number of aminoacids in all peptides/number of aminoacids in protein of interest

    Parameters
    ----------
    protein_seq: str, total sequence of protein of interest
    peptides: list, peptides sequences (obtained from digestion which have lengths in the specified range)
    
    Returns
    -------
    float, protein sequence coverage 
    """
    peptide_len = 0
    for pep in peptides:
        peptide_len+=len(pep)
    coverage = peptide_len/len(protein_seq)*100

    return coverage