def predict_lc_retention_times(peptides):
    """
    This function maps peptide sequences to their predicted retention time using the achrom model RCs_guo_ph7_0.
    
    Parameters
    ----------
    peptides: list, a list of peptide sequences

    Returns
    -------
    dict, dictionary mapping each peptide to its predicted retention time (rounded to two decimals)
    """
    from pyteomics import achrom
    retention_dict = {}
    
    for pep in peptides:
        relative_RT = achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0)
        relative_RT_rounded = round(relative_RT,2)
        retention_dict[pep]= relative_RT_rounded

    return retention_dict



def plot_retention_time(retention_times, resolution=30):
    """
    Plots the distribution of retention times in a histogram.

    Parameters
    ----------
    retention_times: dict,    
            maps peptides to retention times
    resolution: int,
            default=30

    Returns
    -------

    """
    import matplotlib.pyplot as plt
    plt.hist(retention_times.values(), bins=resolution)
    plt.xlabel("retention time (min)")                          #eventhough it is unitless here
    plt.ylabel("intensity")

    plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    This function filters a dictionary of peptides with their according retention times and 
    only saves peptides with retention time in a specified ranged.


    Parameters
    ----------
    peptide_rt_map: dict, a dictionary mapping peptide sequences to their retention time (rounded to two decimals)
    lower_ret_time: float, lower boundary of retention times (is included in range)
    upper_ret_time: float, upper boundary of retention times (is included in range)

    Returns
    -------
    dict, dictionary mapping peptide sequences to their retention time (rounded to two decimals), but only the ones with retention in the specified range

    """
    rt_map_time_sel = {}
    for pep, retention_time in peptide_rt_map.items():
        if lower_ret_time <= retention_time <= upper_ret_time:
            rt_map_time_sel[pep] = retention_time
    
    return rt_map_time_sel   