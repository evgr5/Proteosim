from proteosim.protein_digestion import (
       enzyme_cleavage_patterns,
       digest_protein_sequence,
       digest_protein_collection,
       compute_sequence_coverage,
   )

def test_digest_protein_sequence():
    dummy_seq_1 = "THISKAPIBARAROLLSSUPERFASTINTOTHETREEHOUSE"        # expected outcome: THISK APIBAR AR OLLSSUPER FASTINTOTHE REEHOUSE (6 fragments, 5 should be selected (1 too small, 1 exactly 5; so edge case)) 
    dummy_seq_2 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"                    # sequence that shouldn't be cut and has a length of 30 (should be included in filtered list) #edgecase
    dummy_seq_3 = "WWWWWWWWWKPWWWWWWWWWWRPWWWWWW"                     # sequence should not be cut by trypsin as it has a p after the K and R
    dummy_seq_4 = "PPPPPPPPPPKAPPPPPPPPRAAAAAAAA"                     # sequence should be cut 2 times by trypsin (3 fragments of appropriate length)
   
    enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
    }
    cleave_pattern = enzyme_cleavage_patterns['Trypsin']
    
    #Tests for Trypsin
    fil_pep_tryp, pep_nr_tryp = digest_protein_sequence(dummy_seq_1, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert pep_nr_tryp == 5
    assert fil_pep_tryp == ['THISK', 'APIBAR', 'OLLSSUPER', 'FASTINTOTHETR', 'EEHOUSE']
    fil_pep_tryp, pep_nr_tryp = digest_protein_sequence(dummy_seq_2, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert fil_pep_tryp == ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"]
    fil_pep_tryp, pep_nr_tryp = digest_protein_sequence(dummy_seq_2, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=29)
    assert pep_nr_tryp == 0
    fil_pep_tryp, pep_nr_tryp = digest_protein_sequence(dummy_seq_3, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert fil_pep_tryp == ["WWWWWWWWWKPWWWWWWWWWWRPWWWWWW"]
    fil_pep_tryp, pep_nr_tryp = digest_protein_sequence(dummy_seq_4, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert fil_pep_tryp == ["PPPPPPPPPPK", "APPPPPPPPR", "AAAAAAAA"]
    
    #Test other enzyme patterns
    #LysC
    cleave_pattern = enzyme_cleavage_patterns['LysC']
    fil_pep_lys_c, pep_nr_lys_c = digest_protein_sequence(dummy_seq_3, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert fil_pep_lys_c == ["WWWWWWWWWK", "PWWWWWWWWWWRPWWWWWW"]
    assert pep_nr_lys_c == 2
    #LysN
    cleave_pattern = enzyme_cleavage_patterns['LysN']
    fil_pep_lys_n, pep_nr_lys_n = digest_protein_sequence(dummy_seq_3, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert fil_pep_lys_n == ["WWWWWWWWW", "KPWWWWWWWWWWRPWWWWWW"]
    assert pep_nr_lys_n == 2
    #ArgC
    cleave_pattern = enzyme_cleavage_patterns['ArgC']
    fil_pep_lys_n, pep_nr_lys_n = digest_protein_sequence(dummy_seq_3, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert fil_pep_lys_n == ["WWWWWWWWWKPWWWWWWWWWWR", "PWWWWWW"]
    assert pep_nr_lys_n == 2

test_digest_protein_sequence()


def test_digest_protein_collection():
    dummy_proteins = {
        "dummy_seq_1" : "THISKAPIBARAROLLSSUPERFASTINTOTHETREEHOUSE",       # expected outcome: THISK APIBAR AR OLLSSUPER FASTINTOTHE REEHOUSE (6 fragments, 5 should be selected (1 too small, 1 exactly 5; so edge case))
        "dummy_seq_2" : "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",                   # sequence that shouldn't be cut and has a length of 30 (should be included in filtered list) #edgecase
        "dummy_seq_3" : "WWWWWWWWWKPWWWWWWWWWWRPWWWWWW",                    # sequence should not be cut by trypsin as it has a p after the K and R
        "dummy_seq_4" : "PPPPPPPPPPKAPPPPPPPPRAAAAAAAA"                     # sequence should be cut 2 times by trypsin (3 fragments of appropriate length)
    }
    enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
    }
    cleave_pattern = enzyme_cleavage_patterns['Trypsin']
    #Tests for Trypsin
    pep_dict_tryp = digest_protein_collection(dummy_proteins, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert pep_dict_tryp["dummy_seq_1"] == ['THISK', 'APIBAR', 'OLLSSUPER', 'FASTINTOTHETR', 'EEHOUSE']
    assert pep_dict_tryp["dummy_seq_2"] == ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"]
    assert pep_dict_tryp["dummy_seq_3"] == ["WWWWWWWWWKPWWWWWWWWWWRPWWWWWW"]
    assert pep_dict_tryp["dummy_seq_4"] == ["PPPPPPPPPPK", "APPPPPPPPR", "AAAAAAAA"]
    # length test
    pep_dict_tryp = digest_protein_collection(dummy_proteins, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=29)
    assert len(pep_dict_tryp["dummy_seq_2"]) == 0

    #Test other enzyme patterns
    #LysC
    cleave_pattern = enzyme_cleavage_patterns['LysC']
    pep_dict_lysc = digest_protein_collection(dummy_proteins, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert pep_dict_lysc["dummy_seq_3"] == ["WWWWWWWWWK", "PWWWWWWWWWWRPWWWWWW"]

    #LysN
    cleave_pattern = enzyme_cleavage_patterns['LysN']
    pep_dict_lysn = digest_protein_collection(dummy_proteins, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert pep_dict_lysn["dummy_seq_3"] == ["WWWWWWWWW", "KPWWWWWWWWWWRPWWWWWW"]
    
    #ArgC
    cleave_pattern = enzyme_cleavage_patterns['ArgC']
    pep_dict_argc = digest_protein_collection(dummy_proteins, cleave_pattern=cleave_pattern, min_pep_len=5, max_pep_len=30)
    assert pep_dict_argc["dummy_seq_3"] == ["WWWWWWWWWKPWWWWWWWWWWR", "PWWWWWW"]

test_digest_protein_collection()


def test_compute_sequence_coverage():
    dummy_prot_seq = "AAAAAAAAAAKAAAAAAAAAKAAAA"
    dummy_peps = ["AAAAAAAAAAK", "AAAAAAAAAK"]

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 21/25*100

test_compute_sequence_coverage()