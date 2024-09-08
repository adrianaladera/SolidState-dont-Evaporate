# cd /Users/adrianaladera/Desktop/MIT/research/MOChAs/relaxation_results/originals/3amino/2ISIF
# scp -r aladera@txe1-login.mit.edu:/home/gridsan/aladera/workspace/MOChAs/originals/3amino/2ISIF/SCF .
# cd /Users/adrianaladera/Desktop/MIT/research/MOChAs/relaxation_results/originals/3amino/3ISIF
# scp -r aladera@txe1-login.mit.edu:/home/gridsan/aladera/workspace/MOChAs/originals/3amino/3ISIF/SCF .
# cd /Users/adrianaladera/Desktop/MIT/research/MOChAs/relaxation_results/originals/3aminoH/2ISIF
# scp -r aladera@txe1-login.mit.edu:/home/gridsan/aladera/workspace/MOChAs/originals/3aminoH/2ISIF/SCF .
# cd /Users/adrianaladera/Desktop/MIT/research/MOChAs/relaxation_results/originals/3aminoH/3ISIF
# scp -r aladera@txe1-login.mit.edu:/home/gridsan/aladera/workspace/MOChAs/originals/3aminoH/3ISIF/SCF .

import os
import subprocess

def get_folder(source, dest, folder, mocha, relax):
    #cd = f"cd {dest}/{mocha}/{relax}"
    scp = f"scp -r {source}/{mocha}/{relax}/{folder} {dest}/{mocha}/{relax}"
    #subprocess.run(cd, shell=True)
    subprocess.run(scp, shell=True)

def send_KPOINTS(source, dest, mocha, relax, db):
    scp = f"scp {source}/{mocha}_{relax} {dest}/{mocha}/{relax}/{db}"
    subprocess.run(scp, shell=True)

def send_INCARS(source, dest, mocha, relax, db):
    scp = f"scp {source}/{db} {dest}/{mocha}/{relax}/{db}"
    subprocess.run(scp, shell=True)

mochas = ["3amino"]
relax = ["3ISIF", "2ISIF"]
folder = ["DOS", "SCF"]#,  "bandstructure"]
dos_band = ["bandstructure"]
local_dest = "/Users/adrianaladera/Desktop/MIT/research/MOChAs/relaxation_results/originals"
local_source = "/Users/adrianaladera/Desktop/MIT/research/CSE_thesis/VASP_shit"
remote = "aladera@txe1-login.mit.edu:/home/gridsan/aladera/workspace/MOChAs/originals"

for m in mochas:
    for r in relax:
        # for f in folder:
        #     get_folder(remote, local_dest, f, m, r)
        for db in dos_band:
            if "band" in db:
                send_KPOINTS(f"{local_source}/k-point_brillouin_zone", remote, m, r, db)
            send_INCARS(f"{local_source}/INCARS", remote, m, r, db)