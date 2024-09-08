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
path = "/Users/adrianaladera/Desktop/MIT/research/MOChas/original_structures/3-methoxy/POSCAR"

structure = Structure.from_file("{}".format(path))

# for a in range(len(structure)):
    # if a not in [7,92,58,]:
    # if str(structure[a].species) == "O1" and a in [82,83,84,85,90,91,92,93]: # atom to replace
        # structure.replace(i=a,species='N')
        # print(str(structure[a].species))
        # coordz = structure[a].coords
        # neighbors = structure.get_neighbors(site = structure[a], r = 2.8)
        # for n in neighbors:
        #     print(n.species)
        #     if str(n.species) == "Se1":
        #         print("your mom")
        #         curr_vector = get_components(n.coords, structure[a].coords) # current distance in Å
        #         new_vector = replace_distance(curr_vector, -0.6) # desired length in Å
        #         # new_vector = replace_distance(curr_vector, -0.15)
 
        #         # replace atoms with temp species not already in molecule
        #         structure.replace(i=n.index, species='S')
        #         # structure.replace(i=n.index, species='C', coords=replace_coords(new_vector, n.coords), coords_are_cartesian=True)

# # removing atoms that make up ligand
structure.remove_species(species='O')
structure.remove_species(species='H')
# structure.remove_species(species='C')

path = "/Users/adrianaladera/Desktop/MIT/research/rhino_structure_design/mochas/"
structure.to(filename = "{}DEMO_skeleton.cif".format(path)) # write to file

print("file written your majesty >:)")