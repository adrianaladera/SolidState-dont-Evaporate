import numpy as np
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

def angle_between_vectors(u, v):
    dot_product = sum(i*j for i, j in zip(u, v))
    norm_u = math.sqrt(sum(i**2 for i in u))
    norm_v = math.sqrt(sum(i**2 for i in v))
    cos_theta = dot_product / (norm_u * norm_v)
    angle_rad = math.acos(cos_theta)
    angle_deg = math.degrees(angle_rad)

    return angle_deg