from utils.packages import *

def effmass_from_vbm_cbm(data, settings, vbm, cbm, tol=1e-4):
    ''' data - effmass.inputs.DataVasp() object
        settings - effmass.inputs.Settings() object
        vbm - valence band maximum
        cbm - conduction band minimum
        tol - tolerance depth for generating segments. Default is 1e-4
        
        Returns: effmass.extrema segments thingy'''
    bk_list = []
    for i, en_list in enumerate(data.energies):
        for j, en in enumerate(en_list):
            for k in [vbm, cbm]:
                diff = abs(k - en)
                if diff < tol:
                    bk_list.append([i,j])

    segments = extrema.generate_segments(settings,data, bk = bk_list)
    outputs.plot_segments(data,settings,segments)

    return segments


def get_effmass_data(segments, data, settings, selected_indices=None, path=None):
    ''' segments - effmass yourmom
        data - effmass.inputs.DataVasp() object
        settings - effmass.inputs.Settings() object
        selected_indices - indices to visualize and retrieve data for effective mass.
                           Defaults to None.
        path - path to save. Defaults to None.'''
    
    data_frame = []

    if selected_indices is None:
        # add some option for None on selected_indices
        pass
    else:
        seggs = selected_indices

    for s in seggs:
        data_frame.append({"band": segments[s].band_type, 
        "index":s,
        "max curvature": max(segments[s].five_point_leastsq_fit()),
        "finite diff": segments[s].finite_difference_effmass(),
            "five-pt LSQ":segments[s].five_point_leastsq_effmass()})
    
    outputs.plot_segments(data,settings,[segments[i] for i in selected_indices])
    df = pd.DataFrame(data_frame)
    if path is not None:
        df.to_csv(f"{path}/figures-and-data/effective_mass.csv")
        # fix this
        os.rename("/Users/adrianaladera/Desktop/MIT/research/CSE_thesis/notebooks/effmass_1.png", f"{path}/figures-and-data/effmass.png")

    return df