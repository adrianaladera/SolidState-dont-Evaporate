from utils.packages import *

# from pymatgen.core import Structure, Lattice
# from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
# from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.electronic_structure.core import OrbitalType
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np
# import subprocess as sp
# import pandas as pd
# import os, math, re

def sort_elements(poscar):
    '''Sorts elements between inorganic elements (metals, chalcogens)
        and organic elements.
        poscar: POSCAR file with VASP input'''
    inorganics = ['Ag', 'Au', 'Hg', 'Cu', 'Te', 'Se', 'S']
    inorgs, orgs = [], []
    with open(poscar) as f:
        lines = f.readlines()
        cunt = 0
        flag = 0
        for line in lines:
            for i in line.split():
                if str(i) in inorganics and cunt > 1:
                    flag = 1
            if flag:
                break
            cunt += 1
        for i in lines[cunt].split():
            if i in inorganics:
                inorgs.append(i)
            else:
                orgs.append(i)
    return inorgs, orgs

def projections_to_array(projections):
    """Converts projection data to numpy array.
    
    projections: output of BandStructure.get_projection_on_elements()
    output: np.array of size [n_spin, n_bands, n_panes * distances, n_elements]
    """
    array = []
    for spin in projections.keys():
        spin_array = []
        for band in projections[spin]:
            band_array = []
            for d in band:
                band_array.append(np.array(list(d.values())))
            spin_array.append(np.stack(band_array, axis=0))
        array.append(np.stack(spin_array, axis=0))
    return np.stack(array, axis=0)

def get_element_projection_keys(projections):
    for spin in projections.keys():
        for band in projections[spin]:
            for d in band:
                return d.keys()
            
def get_band_data(path, bs, savedata=False):
    '''Writes BSVasprun data to csv file'''
    data = bs.as_dict()

    vbm = data["vbm"]["energy"]
    cbm = data["cbm"]["energy"]
    efermi = data["efermi"]
    bandgap = data["band_gap"]["energy"]
    transition = data["band_gap"]["transition"]

    df = pd.DataFrame({"vbm": vbm, "cbm":cbm, "eFermi":efermi, "bandgap": bandgap, "transition":transition},index=[0])

    if savedata:
        if not os.path.exists(f"{path}/figures-and-data"):
            os.mkdir(f"{path}/figures-and-data")
        df.to_csv(f"{path}/figures-and-data/estruct_energy_data.csv")
    
    return df