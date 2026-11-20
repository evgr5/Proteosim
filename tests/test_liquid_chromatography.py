from proteosim.liquid_chromatography import (
       predict_lc_retention_times,
       select_retention_time_window,
   )

def test_predict_lc_retention_times():
    from pyteomics import achrom
    peptides = ["THISISATESTSEQ","ATESTWWW","AAAAAA","LLLLLL"]
    expected = {}
    for pep in peptides:
        expected[pep] = round(achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0),2)
    print(expected)

    actual = predict_lc_retention_times(peptides)
    print(actual)
    assert actual == expected

test_predict_lc_retention_times()

def test_select_retention_time_window():
    peptide_rt_map = {'THISISATESTSEQ': 17.3, 'ATESTWWW': 29.5, 'AAAAAA': 13.2, 'LLLLLL': 54.0}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=15, upper_ret_time=30)       #should return 2 values
    assert selected == {'THISISATESTSEQ': 17.3, 'ATESTWWW': 29.5}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=13.2, upper_ret_time=29.5)   #should return 3 values (as boundaries included in range)
    assert selected == {'THISISATESTSEQ': 17.3, 'ATESTWWW': 29.5, 'AAAAAA': 13.2}
test_select_retention_time_window()