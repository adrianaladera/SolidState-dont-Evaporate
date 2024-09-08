import os
import subprocess as sp

with open("KPOINTS", 'r') as f:
    lines = f.readlines()
    # print(lines[-1:][0].strip().split())
    kp = [int(i)*2 for i in lines[-1:][0].strip().split()]
    f.close()
    os.remove("KPOINTS")
    with open("KPOINTS", 'w') as fuck:
        for l in lines[:-1]:
            fuck.write(l)
        fuck.write("{} {} {}".format(kp[0],kp[1],kp[2]))

os.system("cp /global/homes/a/aladera/workspace/VASP/incars/relax3HF . && mv relax3HF INCAR")
