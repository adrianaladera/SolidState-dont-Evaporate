from pymatgen.core import Lattice, Structure, Site
import numpy as np
import os

def get_components(c0, c1):
    '''Parameters:
        - c0: list of Cartesian coordinates for neighboring atom
        - c1: list of Cartesian coordinates for atom to be replaced 
        Returns:
        -  x, y, z components of vector between both sets of coordinates'''
    return (c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2])

def replace_distance(v, dist): # desired distance in Å
    '''Parameters: 
        - v: vector obtained by get_components()
        - desired distance between atoms in Ångstrom 
        Returns:
        - components of new vector '''
    return (v / np.linalg.norm(v)) * dist

def replace_coords(v, c0):
    '''Parameters:
        - v: new vector obtained with replace_distance()
        - c0: coordinates of neighboring atom
        Returns: 
        - updated coordinates of atom that will replace current atom'''
    return v + c0

# path to MOChA file that you want to modify
path = "/Users/adrianaladera/Desktop/MIT/research/mochas/VASP_calculations/originals/glu3_rtr/2ISIF/H-only-relax/checking_struct/scf/"
mocha = "CONTCAR"

structure = Structure.from_file(f"{path}{mocha}".format(path))

# # removing atoms that make up ligand
structure.remove_species(species='C')
structure.remove_species(species='H')
structure.remove_species(species='O')

# for a in range(len(structure)):
#     if str(structure[a].species) == "F1":
#         # structure.replace(idx=a,species='H')
#         print(str(structure[a].species))
#         coordz = structure[a].coords
#         neighbors = structure.get_neighbors(site = structure[a], r = 2.8)
#         for n in neighbors:
#             print(n.species)
#             if str(n.species) == "S1":
#                 print("your mom")
#                 curr_vector = get_components(n.coords, structure[a].coords) # current distance in Å
#                 new_vector = replace_distance(curr_vector, 1.3) # desired length in Å
#                 # new_vector = replace_distance(curr_vector, -0.15)
 
#                 # replace atoms with temp species not already in molecule
#                 # structure.replace(idx=n.index, species='S')
#                 structure.replace(idx=a, species='H', coords=replace_coords(new_vector, n.coords), coords_are_cartesian=True)


structure.to(filename = f"{path}{mocha}_inorganic.cif") # write to file

print("file written your majesty >:)")