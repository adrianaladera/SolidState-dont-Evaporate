import os
import subprocess

#potcar_dir = "/home/gridsan/aladera/workspace/POTCARz/"
potcar_dir = "/Users/adrianaladera/Desktop/MIT/research/POTCARS_PBE.54/"
#poscar_dir = input("Type the dir for the POSCAR: ")
poscar_dir = os.getcwd()
file = input("Input file: ")

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
