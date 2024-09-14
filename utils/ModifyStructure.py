from pymatgen.core import Lattice, Structure, Site, IStructure
from quickmaffs import Vectors
import numpy as np
import os
import ase.geometry
from pymatgen.io.vasp import Poscar
from ase.geometry.analysis import Analysis
import numpy as np
import os
from ase.io import read, write

def replace_ligand_with_H(dest, mocha):
    '''Replaces the ligand at the chalcogen site with a H. 
        Purpose is to isolate the inorganic interaction of the 
        MOCha at the band gap. 
        dest - where to store final file
        mocha - Pymatgen Structure object'''

    # distance in Angstrom between the given chalcogen and a hydrogen
    chalcogen_hydrogen_distances = {'Se': 1.46, 'S': 1.3, 'Te': 1.68}
    chalcogen_carbon_dstances = {'Se': 1.96, 'S': 1.83, 'Te': 2.16}

    all_inorganics = ["Cu", "Au", "Ag", "S", "Te", "Se"]

    species = set([str(s.symbol) for s in mocha.species])

    inorganics = list(set(all_inorganics).intersection(species))
    organics = list(set(species) - set(inorganics))

    for atom in mocha:
        if str(atom.specie) in chalcogen_hydrogen_distances.keys():
            neighbors = mocha.get_neighbors(site = atom, r = chalcogen_carbon_dstances[str(atom.specie)]+0.25)
            if len(neighbors) > 1:
                print("Warning: There is more than one C atom found neighboring your chalcogen. Be sure that there is only one C atom you're looking for!")
            curr_vector = Vectors.get_components(atom.coords, neighbors[0].coords) # current distance in Å
            new_vector = Vectors.replace_distance(curr_vector, chalcogen_hydrogen_distances[str(atom.specie)]) # desired length in Å
    
            # replace atoms with temp species not already in molecule
            mocha.replace(idx=neighbors[0].index, species='Po', coords=Vectors.replace_coords(new_vector, atom.coords), coords_are_cartesian=True)

    #removing atoms that make up ligand
    for o in organics:
        mocha.remove_species(species=o)

    for i, atom in enumerate(mocha):
        if str(atom.specie) == "Po":
            mocha.replace(idx=i, species='H')

    mocha.to(filename = f"{dest}/inorganic-with-H.cif") # write to file
    struct = Structure.from_file(f"{dest}/inorganic-with-H.cif")
    poscar = Poscar(struct)
    poscar.write_file(f"{dest}/inorganic-with-H.vasp")
    os.system(f"rm {dest}/inorganic-with-H.cif") 

def remove_dupes(structure):
    '''Removes duplicates from a structure
        - structure: ase Atoms object'''
    
    print(f"Total: {len(structure)}")
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=True)
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=False)
    print(f"Total with duplicates removed: {len(new_struct)}")

    return structure 

def selective_dynamics(structure, species):
    '''Adds selective dynamics to desired species
        - structure: pymatgen.core.structure object
        - savepath: str filename to save Structure object
        - species: list of strings; species to add selective dynamics'''
    
    select_dicks = []

    for ass in range(len(structure)):
        if str(structure[ass].specie) in species:
            select_dicks.append([True, True, True])
        else:
            select_dicks.append([False, False, False])
    # convert to poscar
    pos = Poscar(structure, selective_dynamics=select_dicks)

    return pos

def remove_species(structure, species):
    '''saves just the inorganic components of a structure to a file
        - mocha: pymatgen.core.Structure object
        - species: string list of species to remove from Structure'''
    
    for s in species:
        structure.remove_species(species=s)
    
    return structure

def find_unit_cell(supercell):
    '''If you have a supercell, get the basic primitive structure
        (unit cell)
        
        - supercell: pymatgen.core.Structure'''
    
    unit_cell = supercell.get_primitive_structure()
    return unit_cell

def modify_lattice(path, dest, file, lattice):
    '''Updates the unit cell of the current .xyz file and saves in .vasp format
    
        - path: the path to the POSCAR
        - dest: where the final converted file will be saved
        - file: name of the .xyz file to be converted
        - lattice: if a 1D array, then lattice will be constructed from parameters
                    a, b, c, alpha, beta, gamma. If 2D, then lattice will be constructed
                    from the given 2D matrix.'''
    if os.path.exists(f'{path}{file}'):
        old_struct = read(f'{path}{file}') 
        structure = remove_dupes(old_struct)
        
        if len(np.array(lattice).shape) == 1: # from 1D array
            lattice = Lattice.from_parameters(lattice[0],lattice[1],lattice[2],lattice[3],lattice[4],lattice[5]) 
            matrix = np.array(lattice)

            new_lattice = Lattice(np.array([matrix[0] - matrix[3],
                                            matrix[1] - matrix[3],
                                            matrix[2] - matrix[3]]))
            structure.set_cell(new_lattice.matrix)  

        else: # from 2D array
            structure.set_cell(lattice) 
        
        write(f"{dest}{file[:-4]}.vasp", structure)
        print(f"Updating unit cell of {file} and writing to *.vasp\n")
        print(f"Written to:\n\n{dest}{file[:-4]}.vasp")