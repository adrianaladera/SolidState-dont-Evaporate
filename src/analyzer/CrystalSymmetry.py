from utils.packages import *

def crystal_system_analyzer(structure):
    ''' Checks to see if the crystal system found
        via Pymatgen's SpacegroupAnalyzer matches
        the definition for the seven crystal systems.
        Returns True if the Pymatgen system and the
        definition system match, False if not.

        structure - pymatgen.core.Structure object'''
    spganal = SpacegroupAnalyzer(structure)
    pmg_crystal = spganal.get_crystal_system()

    a,b,c = structure.lattice.abc
    alpha = structure.lattice.alpha
    beta = structure.lattice.beta
    gamma = structure.lattice.gamma

    if a == b == c and alpha == beta == gamma == 90:
        self_crystal = "cubic"
    elif a == b and a != c and alpha == beta == gamma == 90:
        self_crystal = "tetragonal"
    elif a != b != c and alpha == beta == gamma == 90:
        self_crystal = "orthorhombic"
    elif a == b == c and alpha == beta == gamma != 90:
        self_crystal = "rhombohedral"
    elif a == b != c and alpha == beta == 90 and gamma == 120:
        self_crystal = "hexagonal"
    elif a != b != c and alpha == gamma == 90 != beta:
        self_crystal = "monoclinic"
    elif a != b != c and alpha != beta != gamma:
        self_crystal = "triclinic"
    
    if pmg_crystal != self_crystal:
        print("WARNING: Pymatgen did not determine the correct crystal system. \nYou may need to input the spacegroup manually.")
        print("Lattice parameters: ")
        print(a, b, c, alpha, beta, gamma)
        match = False
    else:
        print(pmg_crystal)
        match = True