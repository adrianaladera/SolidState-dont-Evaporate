from ase.io import read, write
from ase.visualize import view
from ase.visualize.plot import plot_atoms
from ase.io.vasp import write_vasp
import ase.geometry
import os
import subprocess

path = "/Users/adrianaladera/Desktop/MIT/research/mochas/original_structures/1D_Flourescent/AgS_1-naphthyl_1D_Lum/"

def remove_dupes(path, ass):
    np = f"{path}{ass}"
    if os.path.exists(f"{ath}{ass}/"):
        ring = read(np)
        print(f"Total: {len(ring)}")
        dupes = ase.geometry.get_duplicate_atoms(ring, cutoff=0.2, delete=True)
        dupes = ase.geometry.get_duplicate_atoms(ring, cutoff=0.2, delete=False)
        print(f"Total with duplicates removed: {len(ring)}")
        # if len(dupes) == 0:
        #     print("EMPTY! >:D")
        # write_vasp(f"{np}POSCAR/_test", ring)
        write(f"{dest}{ass[:-4]}/test.xyz", ring)
        subprocess.run(f"head {dest}{ass[:-4]}/test.xyz", shell=True)
        # subprocess.run(f"head {dest}{ass[:-4]}/ring.xyz", shell=True)