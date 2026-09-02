import os
import pandas as pd
import math

def total_force(forces):
    '''Total forces in eV / Angstrom on a site (atom)'''
    total_force = math.sqrt(forces[0]**2 + forces[1]**2 + forces[2]**2)
    net_force = forces[0] + forces[1] + forces[2]

    return total_force, net_force

def read_forces(outcar):
    '''from Github user gVallverdu: gVallverdu/read_forces.py'''
    forces = outcar.read_table_pattern(
    header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
    row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
    footer_pattern=r"\s--+",
    postprocess=lambda x: float(x),
    last_one_only=False)

    return forces

def get_forces(path, outcar, structure, inorganics):
    force_info_list = []
    df = pd.DataFrame(columns=["species", "total force (ev/Å)", "net force (ev/Å)", "F_x", "F_y", "F_z"])
    forces = read_forces(outcar)

    # getting the total force eV/Angstrom per atom in the inorganic site
    for atom, f, in zip(range(len(structure)), forces[0]):
        tot_force, net_force = total_force(f)
        for i in inorganics:
            if str(structure[atom].species)[:-1] == i:
                print(f"{str(structure[atom].species)[:-1]} - TOTAL FORCE: {round(tot_force,4)} ev/Å\tNET FORCE: {round(net_force,4)} ev/Å")
                force_info_list.append({"species":str(structure[atom].species)[:-1], 
                           "total force (ev/Å)":tot_force, 
                           "net force (ev/Å)":net_force, 
                           "F_x":f[0], 
                           "F_y":f[1], 
                           "F_z":f[2]})
                
    df = pd.DataFrame(force_info_list)
    if not os.path.exists(f"{path}/figures-and-data"):
        os.mkdir(f"{path}/figures-and-data")
    df.to_csv(f"{path}/figures-and-data/forces_on_inorganics.csv")