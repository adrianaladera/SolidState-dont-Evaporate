from pymatgen.core import Lattice, Structure, Site
import numpy as np
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

path = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/original_pre-relaxed/cube_octahedron/"
if os.path.exists("{}mithreneH_cube_oct_shift.vasp".format(path)):
    os.rename("{}mithreneH_cube_oct_shift.vasp".format(path), "{}POSCAR".format(path))
    print("your mom")

structure = Structure.from_file("{}POSCAR".format(path))

fig = plt.figure()
ax = plt.figure().add_subplot(projection='3d')

for s in structure:
    c = s.coords
    if str(s.species) == "Ag1":
        ax.scatter(c[0], c[1], c[2], color='black')
    if str(s.species) == "Se1":
        ax.scatter(c[0], c[1], c[2], color='green')
    if str(s.species) == "H1":
        ax.scatter(c[0], c[1], c[2], color='pink')

plt.show()