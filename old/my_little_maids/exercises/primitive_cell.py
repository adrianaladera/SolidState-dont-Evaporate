from pymatgen.core import Lattice, Structure, Site, IStructure
import numpy as np
import os
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/original_pre-relaxed/cube_octahedron/"
if os.path.exists(f"{path}mithreneH_cube_oct_c-lattice.vasp"):
    os.rename(f"{path}mithreneH_cube_oct_c-lattice.vasp", f"{path}POSCAR")
    print("your mom")

structure = Structure.from_file("{}POSCAR".format(path))
for s in structure:
    print(s.coords)


lettuce = [[5.938,0,0],[0,7.325,0],[0,0,29.202]]
speecheez = [s.species for s in structure]
coordz = [s.coords for s in structure]
print(len(speecheez), len(coordz))

new_unit_cell = IStructure(lattice=lettuce, species=speecheez, coords=coordz, coords_are_cartesian=True)

new_unit_cell.to(filename = f"{path}test_cube_cell.cif") # write to file
print("file written your majesty >:)")