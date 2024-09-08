import os
import subprocess as sp

for file in os.listdir(os.getcwd()):
    # print(file)
    if "POSCAR-" in file:
        print(file[-3:])
        sp.run("mkdir mi{}".format(file[-3:]), shell=True)
        sp.run("cp {} INCAR POTCAR KPOINTS run.slurm mi{}".format(file, file[-3:]), shell=True)
        sp.run("mv mi{}/{} mi{}/POSCAR".format(file[-3:], file, file[-3:]), shell=True)

sp.run("chmod +x jobcheck.sh", shell=True)
sp.run("./jobcheck.sh", shell=True)

# import subprocess as sp

# def cp_incar():
#     incar_dir = "/global/homes/a/aladera/workspace/VASP/incars/"
#     which_inc = input("\nChoose which incar: \ns - self-consistent\nd - density of states \nb - bandstructure \nr2 - relaxation (2ISIF\n r3 - relaxation (3ISIF)")

#     if which_inc == "s":
#             sp.run(f"cp {incar_dir}scf . && mv scf INCAR", shell = True)                                    
#     elif which_inc == "d":
#             sp.run(f"cp {incar_dir}dos . && mv dos INCAR", shell = True)
#     elif which_inc == "b":
#             sp.run(f"cp {incar_dir}band . && mv band INCAR", shell = True)
#     elif which_inc == "r2":
#             sp.run(f"cp {incar_dir}relax2 . && mv relax2 INCAR", shell = True)
#     elif which_inc == "r3":
#             sp.run(f"cp {incar_dir}relax3 . && mv relax3 INCAR", shell = True)

# def prep_dos(root_dir):
#        '''- root_dir: base of tree , R, where R has ${R}/SCF, ${R}/DOS, ${R}/BAND'''
#        sp.run(f"cp {root_dir}SCF/* {root_dir}DOS" ,shell=True)
#        sp.run(f"rm {root_dir}DOS/INCAR" ,shell=True)
#        cp_incar()


# incar_dir = "/global/homes/a/aladera/workspace/VASP/incars/"