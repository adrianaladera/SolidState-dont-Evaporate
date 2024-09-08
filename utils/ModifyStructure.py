from pymatgen.core import Lattice, Structure, Site, IStructure
from quick_maffs import Vectors
import numpy as np
import os
import ase.geometry
from pymatgen.io.vasp import Poscar
from ase.geometry.analysis import Analysis
import numpy as np
import os
from ase.io import read, write

def just_the_inorganic(structure):
    '''saves just the inorganic components of a structure to a file
        - structure: pymatgen.core.structure object
        - dest: where to save the file to'''
    inorganics_list = ["Au", "Cu", "Ag", "S", "Te", "Se"]
    rm_list = []

    for a in range(len(structure)):
        if str(structure[a].species)[:-1] not in inorganics_list and str(structure[a].species)[:-1] not in rm_list:
            rm_list.append(str(structure[a].species)[:-1])
    
    for r in rm_list:
        structure.remove_species(species=r)

    return structure

def replace_ligand_with_H(dest, structure, removals, temp_site):
    '''Replaces the ligand at the chalcogen site with a H. 
        Purpose is to isolate the inorganic interaction of the 
        structure at the band gap. 
        dest - where to store final file
        structure - Pymatgen Structure object
        removals - list of strings in which each string are the species
        that make up the ligand
        temp_site - site that will be replaced with a hydrogen at the
        appropriate distance from the calchogen'''

    chalcogens = ['Se', 'Te', 'S']

    # distance in Angstrom between the given chalcogen and a hydrogen
    chalcogen_hydrogen_distances = {'Se': 1.46, 'S': 1.3, 'Te': 1.68}

    # removing atoms that make up ligand
    for r in removals:
        structure.remove_species(species=r)

    for a in range(len(structure)):
        if str(structure[a].species) == f"{temp_site}1":    
            print(str(structure[a].species))
            coordz = structure[a].coords
            neighbors = structure.get_neighbors(site = structure[a], r = 2.8)
            for n in neighbors:
                print(n.species)
                if str(n.species)[:-1] in chalcogens:
                    print("your mom")
                    curr_vector = get_components(n.coords, structure[a].coords) # current distance in Å
                    new_vector = replace_distance(curr_vector, chalcogen_hydrogen_distances[str(n.species)[:-1]]) # desired length in Å
    
                    # replace atoms with temp species not already in molecule
                    structure.replace(idx=a, species='H', coords=replace_coords(new_vector, n.coords), coords_are_cartesian=True)

    structure.to(filename = f"{dest}/inorganic-with-H.cif") # write to file
    print(f"file written to {dest}/inorganic-with-H.cif")   

def remove_dupes(structure):
    '''Removes duplicates from a structure
        - structure: ase Atoms object'''
    
    print(f"Total: {len(structure)}")
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=True)
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=False)
    print(f"Total with duplicates removed: {len(new_struct)}")

    return structure 

def selective_dynamics(yourmom, species):
    '''Adds selective dynamics to desired species
        - yourmom: pymatgen.core.structure object
        - filename: str filename to save Structure object
        - species: list of strings; species to add selective dynamics'''
    
    select_dicks = []

    for ass in range(len(yourmom)):
        if str(yourmom[ass].species)[:-1] in species:
            select_dicks.append([True, True, True])
        else:
            select_dicks.append([False, False, False])
    # convert to poscar
    pos = Poscar(yourmom, selective_dynamics=select_dicks)

    return pos

def remove_species(structure, species):
    '''saves just the inorganic components of a structure to a file
        - structure: pymatgen.core.structure object
        - dest: where to save the file to
        - species: string list of species to remove from Structure'''
    
    for s in species:
        structure.remove_species(species=s)

    return structure
    

def find_unit_cell(supercell):
    unit_cell = supercell.get_primitive_structure()
    return unit_cell