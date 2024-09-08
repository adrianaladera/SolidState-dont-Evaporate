import os
import subprocess
import sys

potcar_dir = "/global/homes/a/aladera/workspace/INOLONGERSMOKEPOT/potpaw_PBE/"
# potcar_dir = "/Users/adrianaladera/Desktop/MIT/research/POTCARS_PBE.54/" #local
poscar_dir = os.getcwd()
file = sys.argv[1]

with open("{}/{}".format(poscar_dir, file)) as f:
    print("{}/{}".format(poscar_dir, file))
    lines = f.readlines()
    command = "cat "
    for atom in lines[5].split():
        print(atom)
        command += "{}{}/POTCAR ".format(potcar_dir, atom)
    command += "> {}/POTCAR".format(poscar_dir)
    print(command)
    subprocess.call(command, shell=True)