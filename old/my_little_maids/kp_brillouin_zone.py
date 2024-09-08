from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathLatimerMunro, KPathBase, KPathSetyawanCurtarolo, KPathSeek
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Poscar
import os

def get_kpoint_brillouin_zone(path, dest, mocha, rel_type):
    struct = Structure.from_file(f"{path}scf/CONTCAR")
    # struct = Structure.from_file(f"{path}POSCAR")
    spg_analy =SpacegroupAnalyzer(struct)
    prim_struct=spg_analy.get_primitive_standard_structure(international_monoclinic=False) # based off of Curtarolo convention
    pos = Poscar(prim_struct)
    if not os.path.exists(f"{path}/band_curtarolo"):
        os.mkdir(f"{path}/band_curtarolo")
    # if not os.path.exists(f"{path}/scf_band"):
    #     os.mkdir(f"{path}/scf_band")
    # pos.write_file(f"{path}band/POSCAR_PRIMITIVE.vasp")
    # pos.write_file(f"{path}scf_band/POSCAR_PRIMITIVE.vasp")

    x = input("\nChoose high-symmetry k-path convention (methods can be found on https://pymatgen.org/pymatgen.symmetry.html#module-pymatgen.symmetry.kpath):\n\nl: Latimer-Munro\ns: Setyawan-Curtarolo\nh: Hinuma\n\nYour option: ")
    if x == 'l':
        kpath = KPathLatimerMunro(prim_struct)
        meth = "LatiMunro"
    elif x == 's':
        # kpath = KPathSetyawanCurtarolo(prim_struct)
        kpath = HighSymmKpath(prim_struct) #this is the same as the setyawan thing
        meth = "SetyaCurta"
    elif x == 'h':
        kpath = KPathSeek(prim_struct)
        meth = "Hinuma"
    # kpath = HighSymmKpath(prim_struct)
    else:
        print("Heathen! You did not choose a method! I will choose Latimer-Munro for you.\n\n")
        kpath = KPathLatimerMunro(prim_struct)

    kpts = Kpoints.automatic_linemode(divisions=10,ibz=kpath)
    kpts.write_file(f"{dest}k-point_brillouin_zone/{mocha}_{rel_type}_{meth}")
    kpts.write_file(f"{path}band_curtarolo/KPOINTS")
    print(f"KPOINTS written to {path}band_curtarolo/KPOINTS")
    print(f"Primitive structure written to {path}band_curtarolo/POSCAR_PRIMITIVE.vasp")

# root = f"/Users/adrianaladera/Desktop/MIT/research/mochas/VASP_calculations/originals/1D_Flourescent/{mocha}/"
root = "/Users/adrianaladera/Desktop/MIT/research/mochas/VASP_calculations/fluorescent_1D/Cu_2MMB_ORANGE/"
destination = "/Users/adrianaladera/Desktop/MIT/research/CSE_thesis/VASP_shit/"
stuff = root.split('/')
mocha = stuff[-3] # MUST HAVE '/' AT THE END OF ROOT
relax = stuff[-2]

get_kpoint_brillouin_zone(root, destination, mocha, relax)