from pymatgen.core import Lattice, Structure, Site
import numpy as np
import os
import math

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

def calculate_final_coordinates(initial_coords, length, direction):
    # initial_coords is a tuple (x1, y1, z1)
    x1, y1, z1 = initial_coords

    # direction is a tuple (theta, phi) specifying angles in radians
    theta, phi = direction

    # Calculate the components of the vector
    x_component = length * math.sin(phi) * math.cos(theta)
    y_component = length * math.sin(phi) * math.sin(theta)
    z_component = length * math.cos(phi)

    # Calculate the final coordinates
    x_final = x1 + x_component
    y_final = y1 + y_component
    z_final = z1 + z_component

    return x_final, y_final, z_final

import math

def calculate_direction(initial_coords, final_coords):
    # initial_coords and final_coords are tuples (x1, y1, z1) and (x2, y2, z2)
    x1, y1, z1 = initial_coords
    x2, y2, z2 = final_coords

    # Calculate the differences in coordinates
    delta_x = x2 - x1
    delta_y = y2 - y1
    delta_z = z2 - z1

    # Calculate the magnitude of the vector
    magnitude = math.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

    # Calculate the angles theta and phi (in radians)
    theta = math.atan2(delta_y, delta_x)
    phi = math.acos(delta_z / magnitude)

    return theta, phi

# Example usage:
initial_coords = (1, 2, 3)
final_coords = (4, 5, 6)

direction = calculate_direction(initial_coords, final_coords)
print("Direction (Theta, Phi):", direction)


# Example usage:
initial_coords = (1, 2, 3)
length = 5.0
direction = (math.pi/3, math.pi/4)

final_coords = calculate_final_coordinates(initial_coords, length, direction)
print("Final Coordinates:", final_coords)


# path to MOChA file that you want to modify
path = "/Users/adrianaladera/Desktop/MIT/research/MOChas/original_structures/glu3_rtr/POSCAR"

structure = Structure.from_file("{}".format(path))

for a in range(len(structure)):
    # print(structure[a].index)
    if str(structure[a].species) == "O1" and a not in range(16,24): 
        print(a, structure[a].species)
#         coordz = structure[a].coords
#         neighbors = structure.get_neighbors(site = structure[a], r = 2.8)
#         for n in neighbors:
#             print(n.species)
#             if str(n.species) == "Se1":
#                 print("your mom")
#                 curr_vector = get_components(n.coords, structure[a].coords) # current distance in Å
#                 new_vector = replace_distance(curr_vector, -0.6) # desired length in Å
#                 # new_vector = replace_distance(curr_vector, -0.15)
 
#                 # replace atoms with temp species not already in molecule
#                 structure.replace(i=n.index, species='S')
#                 # structure.replace(i=n.index, species='C', coords=replace_coords(new_vector, n.coords), coords_are_cartesian=True)

# # # removing atoms that make up ligand
# # structure.remove_species(species='Ag')
# # structure.remove_species(species='S')
# # structure.remove_species(species='')

# structure.to(filename = "{}mithrene_S_amine.cif".format(path)) # write to file

# print("file written your majesty >:)")