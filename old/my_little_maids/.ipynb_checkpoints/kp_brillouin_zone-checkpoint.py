from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathLatimerMunro
from pymatgen.symmetry.bandstructure import HighSymmKpath
import os

def get_kpoint_brillouin_zone(path, dest, mocha, rel_type):
    struct = Structure.from_file(f"{path}glu3_rtr_with_H_FINAL.vasp")

    spg_analy =SpacegroupAnalyzer(struct)
    prim_struct=spg_analy.get_primitive_standard_structure(international_monoclinic=False)

    # kpath = KPathLatimerMunro(prim_struct) # for already primitive structures

    kpath = HighSymmKpath(prim_struct)
    print("paths in first brillouin zone :")
    for p in kpath.kpath["path"]:
        print(p)
    kpts = Kpoints.automatic_linemode(divisions=10,ibz=kpath)
    if not os.path.exists(f"{path}/band"):
        os.mkdir(f"{path}/band")
    # kpts.write_file(f"{dest}k-point_brillouin_zone/GLU_TEST")
    kpts.write_file(f"{dest}k-point_brillouin_zone/{mocha}_{rel_type}")

root = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/original_structures/glu3_rtr/"
# destination = root
destination = "/Users/adrianaladera/Desktop/MIT/research/CSE_thesis/VASP_shit/"

# mochas = ["glu3_rtr", "galac_rtr"]
mochas = ["glu3_rtr_UNRELAXED"]
relax = ["3ISIF"]

for m in mochas:
    # for r in relax:
    # path = f"{root}{m}/scf/"
    path = f"{root}"
    get_kpoint_brillouin_zone(path, destination, m, relax[0])
        # path = f"{root}{m}/{r}"
        # get_kpoint_brillouin_zone(path, destination, m, r)