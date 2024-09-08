import numpy as np
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathLatimerMunro, KPathBase, KPathSetyawanCurtarolo, KPathSeek
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Poscar
import os

def make_kpoints(path, filename):
    with open("{}/{}".format(path, filename), 'r') as f:

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

        with open("{}/KPOINTS".format(path), 'w') as kp:
            kp.write("Regular k-point mesh\n")
            kp.write("0\n")
            kp.write("Gamma\n")
            kp.write("{} {} {}\n".format(N1, N2, N3))
        kp.close()
    f.close()

def get_kpoint_brillouin_zone(path, struct, save=False):
    spg_analy =SpacegroupAnalyzer(struct)
    prim_struct=spg_analy.get_primitive_standard_structure(international_monoclinic=False) # based off of Curtarolo convention
    pos = Poscar(prim_struct)
    if save:
        if not os.path.exists(f"{path}/band"):
            os.mkdir(f"{path}/band")
        if not os.path.exists(f"{path}/scf_band"):
            os.mkdir(f"{path}/scf_band")
        pos.write_file(f"{path}band/POSCAR_PRIMITIVE.vasp")
        pos.write_file(f"{path}scf_band/POSCAR_PRIMITIVE.vasp")

    x = input("\nChoose high-symmetry k-path convention (methods can be found on https://pymatgen.org/pymatgen.symmetry.html#module-pymatgen.symmetry.kpath):\n\nl: Latimer-Munro\ns: Setyawan-Curtarolo\nh: Hinuma\n\nYour option: ")
    if x == 'l':
        kpath = KPathLatimerMunro(prim_struct)
    elif x == 's':
        # kpath = KPathSetyawanCurtarolo(prim_struct)
        kpath = HighSymmKpath(prim_struct) #this is the same as the setyawan thing
    elif x == 'h':
        kpath = KPathSeek(prim_struct)
    # kpath = HighSymmKpath(prim_struct)
    else:
        print("Heathen! You did not choose a method! I will choose Latimer-Munro for you.\n\n")
        kpath = KPathLatimerMunro(prim_struct)

    if kpath is not None:
        kpts = Kpoints.automatic_linemode(divisions=10,ibz=kpath)
        kpts.write_file(f"{path}/band/KPOINTS")
        print(f"KPOINTS written to {path}/band/KPOINTS")
        print(f"Primitive structure written to {path}/band/POSCAR_PRIMITIVE.vasp")
    else:
        print("WARNING: your kpath == None. Select a different kpath method or choose the kpath manually.")


