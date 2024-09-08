import os

os.run("grep \"ENCUT\" scf/OUTCAR > scratch-file")
os.run("grep \"NBANDS\" scf/OUTCAR >> scratch-file")

with open("scratch-file", 'r') as f:
    lines = f.readlines()
    with open("band/INCAR", 'a') as fuck:
        fuck.write("NBANDS = {}\n".format(lines[2].split()[14]))
        fuck.write("ENCUT = {}\n".format(lines[0].split()[2]))
    fuck.close()
    with open("dos/INCAR", 'a') as fuck:
        fuck.write("ENCUT = {}\n".format(lines[0].split()[2]))
    fuck.close()
    
f.close()