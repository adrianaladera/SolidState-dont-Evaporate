import os
import subprocess
import numpy as np
from pymatgen.io.vasp import Kpoints
from pymatgen.core.structure import Structure
import os

def make_kpoints(path):
    with open(path, 'r') as f:

        file = f.readlines()
        a1 = np.array([float(file[2].split()[0]), float(file[2].split()[1]), float(file[2].split()[2])])
        a2 = np.array([float(file[3].split()[0]), float(file[3].split()[1]), float(file[3].split()[2])])
        a3 = np.array([float(file[4].split()[0]), float(file[4].split()[1]), float(file[4].split()[2])])

        V = np.linalg.norm(np.dot(a1, np.cross(a2, a3)))
        len_b1 = np.linalg.norm(((2*np.pi) * np.cross(a2, a3)) / V)
        len_b2 = np.linalg.norm(((2*np.pi) * np.cross(a3, a1)) / V)
        len_b3 = np.linalg.norm(((2*np.pi) * np.cross(a1, a2)) / V)

        N1 = round(len_b1 * (1 / min([len_b1, len_b2, len_b3])))
        N2 = round(len_b2 * (1 / min([len_b1, len_b2, len_b3])))
        N3 = round(len_b3 * (1 / min([len_b1, len_b2, len_b3])))

        with open("KPOINTS", 'w') as kp:
            kp.write("Regular k-point mesh\n")
            kp.write("0\n")
            kp.write("Gamma\n")
            kp.write("{} {} {}\n".format(N1, N2, N3))
        kp.close()
    f.close()

file = input("Enter name of file to make KPOINTS and POTCAR from: ")
make_kpoints(file)

#potcar_dir = "/home/gridsan/aladera/workspace/POTCARz/"
potcar_dir = "/Users/adrianaladera/Desktop/MIT/research/POTCARS_PBE.54/"
#poscar_dir = input("Type the dir for the POSCAR: ")
poscar_dir = os.getcwd()

with open("{}/{}".format(poscar_dir, file)) as f:
    lines = f.readlines()

    command = "cat "
    # print(lines[7].split()[2])
    for atom in lines[5].split():
        print(atom)
        command += "{}{}/POTCAR ".format(potcar_dir, atom)
    command += "> {}/POTCAR".format(poscar_dir)
    print(command)
    subprocess.call(command, shell=True)