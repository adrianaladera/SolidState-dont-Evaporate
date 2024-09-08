import pandas as pd
import os
from pymatgen.core import Structure

path = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/prototype_structures/3-methoxy_benzamine/"

# path = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/prototype_structures/mithrene_S_amine/"
if os.path.exists("{}POSCAR.vasp".format(path)):
    os.rename("{}POSCAR.vasp".format(path), "{}POSCAR".format(path))
    print("changing to POSCAR")

structure = Structure.from_file("{}POSCAR".format(path))
excluded_species = ["Ag1", "F1", "Cl1"]
species, x,y,z = [],[],[],[]
for a in range(len(structure)):
    if str(structure[a].species) not in excluded_species:  # atom to replace
        species.append(str(structure[a].species)[:-1])
        coordz = structure[a].coords
        x.append(coordz[0])
        y.append(coordz[1])
        z.append(coordz[2])

df = pd.DataFrame({'species': species, 'x':x, 'y':y, 'z':z})
df.to_csv(f'{path}POSCAR.csv', index=False)
