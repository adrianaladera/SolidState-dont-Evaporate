from pymatgen.core import Lattice, Structure, Site
import numpy as np
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def make_plane(p1, p2, p3):
    v1 = np.array(p2) - np.array(p1)
    v2 = np.array(p3) - np.array(p1)
    equation = np.cross(v1, v2)
    # equation = np.cross(v1, v2) / np.gcd.reduce(np.cross(v1, v2))
    equation = np.append(equation, -np.dot(p1, equation))
    return equation

def mirror_point(a, b, c, d, x1, y1, z1):
    # r = -2*(a*x1 + b*y1 + c*z1 + d)/(a^2 + b^2 + c^2)
    # xp = x1 + a*r
    # yp = y1 + b*r
    # zp = z1 + c*r

    L = a*a + b*b + c*c
    R = d - (a*x1 + b*y1 +c*z1)
    t = R/L
    xp = x1 + 2*t*a
    yp = y1 + 2*t*b
    zp = z1 + 2*t*c

    # k =(-a * x1-b * y1-c * z1-d)/float((a * a + b * b + c * c))
    # x2 = a * k + x1
    # y2 = b * k + y1
    # z2 = c * k + z1
    # xp = 2 * x2-x1
    # yp = 2 * y2-y1
    # zp = 2 * z2-z1
    return [xp, yp, zp]

point = [3,1,2]
ass = make_plane([2,-3,2], [5,-3,-3], [6,3,6])
# a,b,c,d = ass[0], ass[1], ass[2], ass[3]
a,b,c,d = 1,2,1,1
mirror = mirror_point(a, b, c, d, point[0], point[1], point[2])
print(mirror)

X,Y = np.meshgrid(np.linspace(-20,1,20),np.linspace(-20,1,20))
Z = (d - a*X - b*Y) / c

fig = plt.figure()
ax = plt.figure().add_subplot(projection='3d')

ax.plot_surface(X, Y, Z)
ax.scatter(point[0], point[1], point[2], color='red')
ax.scatter(mirror[0], mirror[1], mirror[2], color='black')

plt.show()

path = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/original_pre-relaxed/mithrene/"
if os.path.exists("{}mithrene_edit_cootie_catcher.vasp".format(path)):
    os.rename("{}mithrene_edit_cootie_catcher.vasp".format(path), "{}POSCAR".format(path))
    print("your mom")

structure = Structure.from_file("{}POSCAR".format(path))


ass = make_plane([7.11295, 4.12764, 8.58160], [4.14395, 0.46514, 8.58160], [1.17495, 4.12764, 8.58160])
a,b,c,d = ass[0], ass[1], ass[2], ass[3]

# fig = plt.figure()
# ax = plt.figure().add_subplot(projection='3d')

# for a in range(len(structure)):
#     if str(structure[a].species) != "Se1" and a != 5 and a != 7: 
#         coordz = structure[a].coords
#         mirror_coordz = mirror_point(a, b, c, d, coordz[0], coordz[1], coordz[2])
#         X,Y = np.meshgrid(np.linspace(-15,10,15),np.linspace(-15,10,15))
#         Z = (d - a*X - b*Y) / c 
#         ax.scatter(coordz[0], coordz[1], coordz[2], color='red')
#         ax.scatter(mirror_coordz[0], mirror_coordz[1], mirror_coordz[2], color='black')


# ax.plot_surface(X, Y, Z)
# plt.show()

        # structure.append(species=structure[a].species, coords=mirror_coordz)

# structure.to(filename = "{}mithrene_edit_cube_octahedron.cif".format(path)) # write to file
# print("file written your majesty >:)")

        # neighbors = structure.get_neighbors(site = structure[a], r = 2.15002)
        # plane = []
        # if "C1" and "H1" not in [n.species for n in neighbors] and len(plane) <= 3:
        #     plane.append(coordz)
        # print(plane)


        # if len(set(neighbors))==1 and str(neighbors[0].species) == "Ag1":
#         for n in neighbors:
#             if str(n.species) == "S1":
#                 curr_vector = get_components(n.coords, structure[a].coords) # current distannce in Å
#                 new_vector = replace_distance(curr_vector, 1.38) # desired length in Å
 
#                 # replace atoms with temp species not already in molecule
#                 structure.replace(i=a, species='Cl', coords=replace_coords(new_vector, n.coords), coords_are_cartesian=True)

# # removing atoms that make up ligand
# structure.remove_species(species='C')
# structure.remove_species(species='H')
# structure.remove_species(species='N')

# for a in range(len(structure)): # replacing temp species with desired H species
#     if str(structure[a].species) == "Cl1":
#         structure.replace(i=a, species='H', coords=None)

# structure.to(filename = "{}testing.cif".format(path)) # write to file

# print("file written your majesty >:)")