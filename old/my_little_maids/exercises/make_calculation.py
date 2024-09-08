import subprocess as sp
import os 

mochas = ["3amino", "mithrene"]
relax = ["2ISIF", "3ISIF"]

cwd = input("Enter cwd for calculation: ")

for m in mochas:
    for r in relax:
        path = ""
        if m == "mithrene":
            path = f"{cwd}/{m}/{r}/enc800_kp_541_ediffg-1e-3/"
        else:
            path = f"{cwd}/{m}/{r}"
        if os.path.exists(path):
            print(path)
        # os.chdir(path)
        # with open(f"{path}/relaxation/job.out") as f:
        #     if "reached required accuracy" in f.readlines():
        #         print("your mom")