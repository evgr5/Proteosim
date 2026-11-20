amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    This function calculartes the molecular mass of peptides based on their sequence using the molecular weights of its amino acids.

    Parameters
    ----------
    peptide_seq: str, peptide sequence
    amino_acid_mass_dict: dict, mapping amino acids to molecular weight in dalton

    Returns
    -------
    dict, mapping peptide to mass of this peptide in dalton
    """

    peptide_mass_map = {}
    mol_mass=0

    for amino_acid in peptide_seq:
        mol_mass += amino_acid_mass_dict[amino_acid]
        peptide_mass_map[peptide_seq]= mol_mass
        #print(amino_acid)

    return peptide_mass_map

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    This function calculates the molecular mass of peptides based on their sequence using the 
    molecular weights of its amino acids.

    Parameters
    ----------
    peptide_seq: str, peptide sequence
    amino_acid_mass_dict: dict, mapping amino acids to molecular weight in dalton

    Returns
    -------
    dict, mapping peptide to mass of this peptide in dalton
    """

    peptide_mass_map = {}
    for peptide_seq in peptides:
        mol_mass=0
        for amino_acid in peptide_seq:
            mol_mass += amino_acid_mass_dict[amino_acid]
        peptide_mass_map[peptide_seq]= mol_mass
        

    return peptide_mass_map


def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Calculates the m/z value = (neutral_mass + charge * proton_mass) of a dictionary with peptides 
    and neurtral molecular masses and creates a dictionary of the peptides and their m/z values.
    
    Parameters
    ----------
    peptide_mass_map:   dict,
                        dictionary mapping peptides to their molecular mass

    charge:     float,
                charge of the peptide

    proton_mass:    float, 
                    mass of a proton, default=1.007

    Returns
    -------
    dict, dictionary mapping peptides to their m/z values
    """
    mz_peptide_map={}

    for peptide in peptide_mass_map:
        mz_value = (peptide_mass_map[peptide] + charge * proton_mass) / charge
        mz_peptide_map[peptide]= mz_value

    return mz_peptide_map


import numpy as np
import matplotlib.pyplot as plt
import random

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Add a short description here.
    MS1 spectrum displays peptide ion intensities (y-axis) versus their m/z values (x-axis)

    Parameters
    ----------
    mz_values: list, m/z values
    random_count_range: 
    Returns
    -------

    """
    intensities =[]
    min, max = random_count_range
    for i in range(len(mz_values)):
        np.random.seed(seed)
        intensities.append(random.randint(min, max+1))

    plt.bar(x=mz_values, height=intensities, width=1)


def fragment_peptide(peptide):
    """
    Creates list of possible peptides when peptide is fragmented in two pieces

    Parameters
    ----------
    peptide:    str, peptide sequence

    Returns
    -------
    list, List of possible peptides (both b_ion-series and y_ion_series)

    """
    combi_frag_list = []
    b=0
    y=0

    for b in range (1, len(peptide)):
        combi_frag_list.append(peptide[:b])
        b += 1

    for y in range (1, len(peptide)):
        combi_frag_list.append(peptide[-y:])
        y += 1

    print(combi_frag_list)

    return combi_frag_list
