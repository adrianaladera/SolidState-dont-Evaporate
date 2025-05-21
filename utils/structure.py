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
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import write
import numpy as np
import pubchempy as pcp

def remove_dupes(structure):
    '''Removes duplicates from a structure
        - structure: ase Atoms object'''
    
    print(f"Total: {len(structure)}")
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=True)
    new_struct = ase.geometry.get_duplicate_atoms(structure, cutoff=0.2, delete=False)
    print(f"Total with duplicates removed: {len(new_struct)}")

    return structure 

class ModifyStructure:
    def __init__(self, path, dest, structure):
        self.path = path
        self.dest = dest
        self.structure = structure

    def replace_ligand_with_H(self, tolerance=0.25):
        '''Replaces the ligand at the chalcogen site with a H. 
            Purpose is to isolate the inorganic interaction of the 
            MOCha at the band gap. 
            dest - where to store final file
            mocha - Pymatgen Structure object'''

        # distance in Angstrom between the given chalcogen and a hydrogen
        chalcogen_hydrogen_distances = {'Se': 1.46, 'S': 1.3, 'Te': 1.68}
        chalcogen_carbon_dstances = {'Se': 1.96, 'S': 1.83, 'Te': 2.16}

        all_inorganics = ["Cu", "Au", "Ag", "S", "Te", "Se"]

        dest = self.dest
        mocha = self.structure

        species = set([str(s.symbol) for s in mocha.species])

        inorganics = list(set(all_inorganics).intersection(species))
        organics = list(set(species) - set(inorganics))

        notempty = False
        for atom in mocha:
            if str(atom.specie) in chalcogen_hydrogen_distances.keys():
                neighbors = mocha.get_neighbors(site = atom, r = chalcogen_carbon_dstances[str(atom.specie)]+tolerance)
                if len(neighbors) > 1:
                    print("Warning: There is more than one C atom found neighboring your chalcogen. Be sure that there is only one C atom you're looking for!")
                if len(neighbors) != 0:
                    notempty = True
                    curr_vector = Vectors.get_components(atom.coords, neighbors[0].coords) # current distance in Å
                    new_vector = Vectors.replace_distance(curr_vector, chalcogen_hydrogen_distances[str(atom.specie)]) # desired length in Å
            
                    # replace atoms with temp species not already in molecule
                    mocha.replace(idx=neighbors[0].index, species='Po', coords=Vectors.replace_coords(new_vector, atom.coords), coords_are_cartesian=True)

        #removing atoms that make up ligand
        if notempty:
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
            print(f"file written to {dest}")
        else:
            print("There are no neighbors! Try increasing the tolerance for chalc-carbon distances.")
            print("Default = 0.25 Å")

    def selective_dynamics(self, species):
        '''Adds selective dynamics to desired species
            - structure: pymatgen.core.structure object
            - savepath: str filename to save Structure object
            - species: list of strings; species to add selective dynamics'''
        
        structure = self.structure
        
        select_dicks = []

        for ass in range(len(structure)):
            if str(structure[ass].specie) in species:
                select_dicks.append([True, True, True])
            else:
                select_dicks.append([False, False, False])
        # convert to poscar
        pos = Poscar(structure, selective_dynamics=select_dicks)

        return pos

    def remove_species(self, species):
        '''saves just the inorganic components of a structure to a file
            - mocha: pymatgen.core.Structure object
            - species: string list of species to remove from Structure'''
        
        structure = self.structure
        
        for s in species:
            structure.remove_species(species=s)
        
        return structure


class ModifyLattice:
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

class Ligands:
    def __init__(self, path, save=False):
        self.path = path
        self.save = save

    def build_ligand(self, iupac_name, lattice_dims):
        root = self.path

        compound = pcp.get_compounds(iupac_name, namespace='name')

        # Step 2: Get the SMILES string
        if compound:
            smiles = compound[0].isomeric_smiles
            print(iupac_name + ": " + smiles)

            # Step 1: Generate RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

            # Step 2: Extract positions and elements
            conf = mol.GetConformer()
            positions = np.array([
                [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                for i in range(mol.GetNumAtoms())
            ])
            symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

            # Step 3: Create ASE Atoms object and define a cubic cell
            atoms = Atoms(symbols=symbols, positions=positions)
            atoms.center(vacuum=0.0)  # ensure it's roughly centered
            atoms.set_cell([lattice_dims, lattice_dims, lattice_dims])
            atoms.set_pbc([True, True, True])  # Periodic BCs if needed

            if self.save:
                if not os.path.exists(f"{root}{iupac_name}"):
                    os.mkdir(f"{root}{iupac_name}")

                write(f"{root}{iupac_name}/POSCAR.vasp", atoms, format="vasp")
        else:
            print(iupac_name + " not found 😢")