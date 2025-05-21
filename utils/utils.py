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

METALS = ['Ag', 'Au', 'Hg', 'Cu']
CHALCS = ['Te', 'Se', "S"]

def sort_elements(poscar):
    '''Sorts elements between inorganic elements (metals, chalcogens)
        and organic elements.
        poscar: POSCAR file with VASP input'''
    inorganics = METALS + CHALCS
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