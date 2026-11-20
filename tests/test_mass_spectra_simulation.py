from proteosim.mass_spectra_simulation import (
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz_collection,
    plot_spectrum,
    fragment_peptide,
    amino_acid_mass_dalton,
)

def test_calculate_mol_mass():
    aa_mass_dict = {
        "A":1, "B":2, "C":3
    }
    test_peptide_1 = "AAAAAA"
    test_peptide_2 = "BC"
    test_peptide_3 = ""
    expected_mass_1 = {"AAAAAA":6.00}
    expected_mass_2 = {"BC":5.00}
    expected_mass_3 = {}

    calculated_mass_1 = calculate_mol_mass(test_peptide_1, aa_mass_dict)
    calculated_mass_2 = calculate_mol_mass(test_peptide_2, aa_mass_dict)
    calculated_mass_3 = calculate_mol_mass(test_peptide_3, aa_mass_dict)
    
    assert calculated_mass_1 == expected_mass_1
    assert calculated_mass_2 == expected_mass_2
    assert calculated_mass_3 == expected_mass_3

test_calculate_mol_mass()


def test_calculate_mol_mass_collection():
    aa_mass_dict = {"A": 1, "P": 2, "E": 3, "W": 4}
    peptides = ["AAA","PEP","W"]
    expected = {"AAA": 3,"PEP": 7,"W":4}
    actual = calculate_mol_mass_collection(peptides, aa_mass_dict)

    assert actual == expected

test_calculate_mol_mass_collection()


def test_calculate_mz_collection():
    peptide_mass_map = {"AAA": 3,"PEP": 7,"W":4}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {"AAA": (3+2.014)/2,"PEP": (7+2.014)/2,"W":(4+2.014)/2}
#(neutral_mass + charge * proton_mass) / charge
#(neutral_mass + 2 * 1.007) / 2      2*1.007=2.014
    assert actual == expected

test_calculate_mz_collection()


def test_fragment_peptide():
    peptide = "PEPT"
    expected = ["P", "PE", "PEP", "T" , "PT", "EPT"]
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)
    assert "PEPT" not in expected  # which string should not be in the expected?

test_fragment_peptide()