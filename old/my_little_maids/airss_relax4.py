import time
import subprocess as sp
import os

def is_converged(path, seed, ionic_steps):
    '''ionic_steps - number of steps (AIRSS needs only 3 or 4 steps, considered "converged")'''
    sp.run("grep \"{} F=\" {}/job.out | awk '{}print $1{}' > {}/ionic_steps".format(ionic_steps, path, '{','}', path), shell=True)
    with open("{}/ionic_steps".format(path), "r") as f:
        if str(ionic_steps) in f.read():
            return 1
        else:
            return 0
        
def is_running(path, dir_name):
    '''Checks to see if a job is currently running in a directory. If so, the job is not overridden'''
    # os.chdir(path)
    os.system("squeue -u aladera --format \"%Z\" > running_jobs")
    with open("running_jobs", 'r') as f:  
        if dir_name in f.read():
            return 1
        else:
            return 0

path = "/pscratch/sd/a/aladera/airss/AgS/var1/"
cunt = 1
while cunt < 6:
    print("Round {}".format(cunt))
    time.sleep(14400)
    converged = []
    unconverged = []
    for i in os.listdir(path):
        if os.path.isdir(i) and seed in i:
            if is_converged("{}{}".format(path, i), seed, 4):
                converged.append(i)
            else:
                unconverged.append(i)       
    print("Converged: {}, Unconverged: {}".format(len(converged), len(unconverged)))
    
    for i in range(len(converged)):
        print(converged[i])
        os.chdir(path + converged[i])
        if is_running(path, converged[i]):
            print("{} already in queue".format(converged[i]))
            pass
        else:
            os.system("cp CONTCAR POSCAR")
            for j in range(1,5):
                if not os.path.exists("POSCAR_{}".format(j)):
                    print(converged[i], j)
                    os.system("cp CONTCAR POSCAR_{}".format(j))
                    break
            os.system("vrun 1 32 5 {}".format(converged[i]))
            os.chdir(path)
        cunt += 1

