from pymatgen.core import Lattice, Structure, Site
import numpy as np
import os
import numpy as np
import matplotlib.pyplot as plt

def tetrahedron_centroid(p1,p2,p3, p4):
    px = (p1[0] + p2[0] + p3[0] + p4[0])/4
    py = (p1[1] + p2[1] + p3[1] + p4[1])/4
    pz = (p1[2] + p2[2] + p3[2] + p4[2])/4

    return [px, py, pz]
        

path = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/original_pre-relaxed/cube_octahedron/"
if os.path.exists("{}stuff_14_your mom.vasp".format(path)):
    os.rename("{}stuff_14_your mom.vasp".format(path), "{}POSCAR".format(path))
    print("your mom")

structure = Structure.from_file("{}POSCAR".format(path))

for i,s in enumerate(structure):
    print(i, s.coords)

tetrahedra = [[2,4,3,5],[2,11,1,5],[9,11,5,8],[4,9,10,5],
              [12,0,3,5],[0,6,1,5],[6,7,8,5],[7,12,10,5]]

Ag_indices = []

for s in range(len(structure)):
    structure.replace(i=s, species='Se')
# print(len(structure))

4:0.10000   0.10000   0.50000
2:0.10000   0.50000   0.90000

# Ag_ind = tetrahedron_centroid(structure[4].coords, structure[9].coords, structure[10].coords, structure[5].coords)
# structure.append(species='Ag', coords=Ag_ind)

# for i,t in enumerate(tetrahedra):
#     print(i)
#     if i == 0:
#         structure = Structure.from_file("{}POSCAR".format(path))
#     else:
#         structure = Structure.from_file("{}cube_octahedron.cif".format(path))
#     Ag_ind = tetrahedron_centroid(structure[t[0]].coords, structure[t[1]].coords, structure[t[2]].coords, structure[t[3]].coords)
#     print(Ag_ind)
#     structure.append(species='Ag', coords=Ag_ind)
#     structure.to(filename = "{}cube_octahedron.cif".format(path)) # write to file

structure.to(filename = "{}cube_octahedron.cif".format(path)) # write to file
print("file written your majesty >:)")