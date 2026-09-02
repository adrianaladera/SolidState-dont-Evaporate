from ase.io import read, write
from ase.build.molecule import molecule
from ase.io import write
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
import os
import math
import numpy as np
from ase import Atoms
from pymatgen.core.lattice import Lattice
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
import subprocess as sp
import ase.geometry

############################### GLOBALS ###############################

# for POSCAR to XYZ-RGB
elements_to_rgb = {'Ag':[192,192,192], 'Hg':[210, 178, 221], 'Au':[255, 179, 0], 'Cu':[6, 9, 189],
                   'Te':[151, 136, 72], 'S':[255,255,51], 'Se':[57, 255, 20],
           'C':[102,0,0], 'N':[102,255,255], 'H':[255,204,204], 'O':[255,0,0], 'F':[130,197,243], 
           'Br':[94, 35, 1],'Na':[209, 185, 2],
           'origin':[255,255,255], 'point':[0,0,0]}

# for XYZ-RGB to POSCAR
rgb_to_elements = {}
for s in elements_to_rgb:
    rgb = elements_to_rgb[s]
    rgb_to_elements[f"{rgb[0]}.{rgb[1]}.{rgb[2]}"] = s


############################### FUNCTIONS ###############################

def poscar_to_xyz_rgb(path, filename):
    '''Converts a VASP file into a .xyz readable points file for Rhino, with elements
        mapped to RGB values
        
        - path: the path to the POSCAR
        - filename: name of the file to be converted(without the extension)
        - extension: .vasp or '' if POSCAR'''
    poscar = Structure.from_file(f"{path}{filename}")
    run = True
    while run:
        if os.path.exists(f"{path}/xyz-rgb-points.xyz"):
            sp.run(f"rm {path}/xyz-rgb-points.xyz", shell=True)
            print("\nFile already exists. Rewriting now...")
            
        with open(f"{path}/xyz-rgb-points.xyz", 'a') as f:
            for i in poscar.as_dict()['sites']:
                colors = species_to_rgb(str(i['species'][0]['element']))
                lines = [f"{i['xyz'][0]} {i['xyz'][1]} {i['xyz'][2]} {colors[0]} {colors[1]} {colors[2]}\n"]
                f.writelines(lines)
            f.close()
            print("\nConverted POSCAR to XYZRGB points file")
            print(f"Written to: \n{path}/xyz-rgb-points.xyz")
            run = False
        
def species_to_rgb(species, rgb=elements_to_rgb):
    return rgb[species]

def remove_dupes(structure):
    print(f"Total: {len(structure)}")
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=True)
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=False)
    print(f"Total with duplicates removed: {len(new_struct)}")

    return structure

def xyz_rgb_to_poscar(path, dest, filename, extension):
    '''Converts a Rhino .txt readable points back into an .xyz file, with RGB values
        mapped to back to elements
        
        - path: the path to .txt file
        - dest: where the final converted file will be saved
        - filename: name of the file to be converted(without the extension)
        - extension: .txt'''
    run = True
    flag = 0
    lattice_points = []
    lattice_origin = None
    while run:
        if not os.path.exists(f"{dest}{filename}_FINAL.xyz") or flag == 1:
            with open(f"{path}{filename}{extension}", 'r') as f:
                with open(f"{dest}{filename}_FINAL.xyz", 'w') as fuck:
                    lines = f.readlines()
                    species_dict, species_list = {},{}
                    for i in lines:
                        i = i.strip('\n')
                        species = rgb_to_species([int(i.split(',')[3]), int(i.split(',')[4]), int(i.split(',')[5])])
                        if species not in species_dict.keys() and species != "origin" and species != "point":
                            species_dict[species] = f"{int(i.split(',')[3])} {int(i.split(',')[4])} {int(i.split(',')[5])}"
                            species_list[species] = []
                    for i in lines:
                        i = i.strip('\n')
                        species = rgb_to_species([int(i.split(',')[3]), int(i.split(',')[4]), int(i.split(',')[5])])
                        if species != "origin" and species != "point":
                            species_list[species].append([float(i.split(',')[0]), float(i.split(',')[1]),float(i.split(',')[2])])
                        elif species == "point":
                            lattice_points.append([float(i.split(',')[0]), float(i.split(',')[1]),float(i.split(',')[2])])
                        elif species == "origin":
                            lattice_origin = [float(i.split(',')[0]), float(i.split(',')[1]),float(i.split(',')[2])]
                    s, s2 = '', ''
                    s += f"{sum([len(species_list[spec]) for spec in species_list])}\n"
                    for spec in species_list:
                        s += f"{spec}{len(species_list[spec])} "
                        for a in species_list[spec]:
                            s2 += f"{spec} {a[0]} {a[1]} {a[2]}\n"
                    s += '\n' + s2
                    fuck.write(s)
            print("Converting text file to XYZ readable molecule file")
            print(f"Written to {dest}{filename}_FINAL.xyz")
            run = False
            fuck.close()
        else:
            sp.run(f"rm {dest}{filename}_FINAL.xyz", shell=True)
            print("\nFile already exists. Rewriting now...")
            run = True
            flag = 1

    lattice_points.append(lattice_origin)
    print(lattice_points)
    return lattice_points
              
def rgb_to_species(rgb_vals, rgb = rgb_to_elements):
    '''Converts rgb values of Rhino points into POSCAR readable species'''
    key = f"{rgb_vals[0]}.{rgb_vals[1]}.{rgb_vals[2]}"
    return rgb[key]