from utils.packages import *

def regular_kpoints(path, structure, factor=None, mesh_type=None):
    '''Create regular k-points based on the lattice vectors given 
        in the POSCAR file.
        path - where KPOINTS is written
        structure - pymatgen.core.Structure object
        factor - scales each corresponding kpoint instance for
                greater kpoint density. Defaults to None.
        mesh_type - Monkhorst-pack or Gamma. Defaults to Gamma.'''
    
    lattice = structure.lattice.matrix
    a1, a2, a3 = lattice[0], lattice[1], lattice[2]

    V = np.linalg.norm(np.dot(a1, np.cross(a2, a3)))
    len_b1 = np.linalg.norm(((2*np.pi) * np.cross(a2, a3)) / V)
    len_b2 = np.linalg.norm(((2*np.pi) * np.cross(a3, a1)) / V)
    len_b3 = np.linalg.norm(((2*np.pi) * np.cross(a1, a2)) / V)

    N1 = round(len_b1 * (1 / min([len_b1, len_b2, len_b3])))
    N2 = round(len_b2 * (1 / min([len_b1, len_b2, len_b3])))
    N3 = round(len_b3 * (1 / min([len_b1, len_b2, len_b3])))

    # multiply basic kpoints by a factor for greater kpoint density
    if factor is not None:
        N1 *= factor
        N2 *= factor
        N3 *= factor

    s = "Regular k-point mesh\n0\n"

    if mesh_type is None: #default
        s += "Gamma\n"
    else:
        s += "Monkhorst-pack\n"
    s += "{} {} {}\n".format(N1, N2, N3)


    with open(f"{path}/KPOINTS", 'w') as kp:
        kp.write(s)
        kp.close()
    
    print("Regular KPOINTS written")

def kpoints_by_density(path, structure, density=1000):
    kpoints = Kpoints.automatic_density(structure, kppa=density)
    kpoints.write_file(f"{path}/KPOINTS")


def band_kpoints(struct, path, kpath_type='s',ishybrid=False):
    ''' structure - pymatgen.core.Structure object
        path - path where the band/ and scf_band/ directories
                will be created
        ishybrid - whether or not the band calculation being
                    created will use a hybrid functional (b3lyp).
                    Defaults to False.'''
    spg_analy =SpacegroupAnalyzer(struct)
    prim_struct=spg_analy.get_primitive_standard_structure(international_monoclinic=False) # based off of Curtarolo convention
    pos = Poscar(prim_struct)
    if not os.path.exists(f"{path}/band"):
        os.mkdir(f"{path}/band")
    if not os.path.exists(f"{path}/scf_band"):
        os.mkdir(f"{path}/scf_band")
    pos.write_file(f"{path}/band/POSCAR")
    pos.write_file(f"{path}/scf_band/POSCAR")

    # x = input("\nChoose high-symmetry k-path convention (methods can be found on https://pymatgen.org/pymatgen.symmetry.html#module-pymatgen.symmetry.kpath):\n\nl: Latimer-Munro\ns: Setyawan-Curtarolo\nh: Hinuma\n\nYour option: ")
    x = kpath_type
    if x == 'l':
        kpath = KPathLatimerMunro(prim_struct)
    elif x == 's':
        kpath = HighSymmKpath(prim_struct) #this is the same as KPathSetyawanCurtarolo(prim_struct)
    elif x == 'h':
        kpath = KPathSeek(prim_struct)
    else:
        print("Heathen! You did not choose a method! I will choose Latimer-Munro for you.\n\n")
        kpath = KPathLatimerMunro(prim_struct)

    if kpath.kpath is not None:
        kpts = Kpoints.automatic_linemode(divisions=6,ibz=kpath)
        if ishybrid:
            kpts.write_file(f"{path}/band/KPOINTS_OPT")
            os.system(f"cp {path}/scf/KPOINTS {path}/band/KPOINTS")
            print("Band KPOINTS and KPOINTS_OPT written for B3LYP functional")
        else:
            kpts.write_file(f"{path}/band/KPOINTS")
            print("Band KPOINTS written for default functional (PBE)")
    
    
    else:
        print("kpath is empty! You may have to determine the spacegroup manually:")
    
def band_incar(potcar_path, write_path, ishybrid=False,
               algo=None, prec=None, npar=None, ncore=None, 
             kpar=None, lreal=None, lcharg=None, lwave=None, ismear=None, sigma=None):
    '''MUST have the outcar path of scf_band
    write_path - must be a directory called scf_band'''
    ass = True
    if ass is True:
        params = {
        "ALGO": algo if algo is not None else "Normal",
        "PREC": prec if prec is not None else "Accurate",
        "NPAR": npar if npar is not None else 8,
        "NCORE": ncore if ncore is not None else 8,
        "KPAR": kpar if kpar is not None else 2,
        "LREAL": lreal if lreal is not None else "Auto",
        "LCHARG": lcharg if lcharg is not None else "False",
        "LWAVE": lwave if lwave is not None else "False",
        "IMSEAR": ismear if ismear is not None else 0,
        "SIGMA": sigma if sigma is not None else 0.03
        }

        params["NELM"] = 250
        params["NELMIN"] = 4
        params["ISTART"] = 1    
        params["EDIFF"] = 1e-6
        params["EDIFFG"] = 1e-4
        params["ISIF"] = 2 
        params["NSW"] = 0
        params["IBRION"] = -1
        params["ISYM"] = 2
        params["ICHARG"] = 11
        params["LORBIT"] = 11
    
        command = f"grep \"ENMAX\" {potcar_path}"
        os.system(command)
        output = os.popen(command).read();
        energy_vals = [float(i) for i in output.replace(";", "").split() if re.match(r'^-?\d+(?:\.\d+)$', i) is not None]
        params["ENCUT"] = 1.3 * max(energy_vals)

        if ishybrid:
            params["LHFCALC"] = ".TRUE." 
            params["GGA"]     = "B3"
            params["AEXX"]    = 0.2
            params["AGGAX"]   = 0.72 
            params["AGGAC"]   = 0.81 
            params["ALDAC"]   = 0.19

        incar = Incar.from_dict(params)
        incar.write_file(f"{write_path}/INCAR")

        print("Band INCAR written")

def dos_incar(potcar_path, write_path, ishybrid=False,
               algo=None, prec=None, npar=None, ncore=None, 
             kpar=None, lcharg=None, lwave=None, lreal=None, ismear=None, sigma=None):
    '''MUST have the outcar path of scf_band
    write_path - must be a directory called scf_band'''
    ass = True
    if ass is True:
        params = {
        "ALGO": algo if algo is not None else "Normal",
        "PREC": prec if prec is not None else "Accurate",
        "NPAR": npar if npar is not None else 8,
        "NCORE": ncore if ncore is not None else 8,
        "KPAR": kpar if kpar is not None else 2,
        "LREAL": lreal if lreal is not None else "Auto",
        "LCHARG": lcharg if lcharg is not None else "False",
        "LWAVE": lwave if lwave is not None else "False",
        "IMSEAR": ismear if ismear is not None else 0,
        "SIGMA": sigma if sigma is not None else 0.03
        }

        params["NELM"] = 250
        params["NELMIN"] = 4
        params["ISTART"] = 1    
        params["EDIFF"] = 1e-6
        params["EDIFFG"] = 1e-4
        params["ISIF"] = 2 
        params["NSW"] = 0
        params["IBRION"] = -1
        params["ISYM"] = 2
        params["ICHARG"] = 11
        params["LORBIT"] = 11
        params["NEDOS"] = 2000

        command = f"grep \"ENMAX\" {potcar_path}"
        os.system(command)
        output = os.popen(command).read();
        energy_vals = [float(i) for i in output.replace(";", "").split() if re.match(r'^-?\d+(?:\.\d+)$', i) is not None]
        params["ENCUT"] = 1.3 * max(energy_vals)

        if ishybrid:
            params["LHFCALC"] = ".TRUE." 
            params["GGA"]     = "B3"
            params["AEXX"]    = 0.2
            params["AGGAX"]   = 0.72 
            params["AGGAC"]   = 0.81 
            params["ALDAC"]   = 0.19

        incar = Incar.from_dict(params)
        if not os.path.exists(f"{write_path}/dos"):
            os.mkdir(f"{write_path}/dos")
            incar.write_file(f"{write_path}/dos/INCAR")

        print("Density of States INCAR written")


def make_potcar(path):
    '''Creates a potcar by concatenating all the pseudoponetials
        (POTCARs) of the species in the structure.
        path - path to save resulting POTCAR. POTCAR will be saved
            to the same directory as the POSCAR.'''
    potcar_dir = "/Users/adrianaladera/Desktop/MIT/research/POTCARS_PBE.54/"

    with open(f"{path}/POSCAR", 'r') as f:
        lines = f.readlines()
        command = "cat "

        for atom in lines[5].split():
            command += "{}{}/POTCAR ".format(potcar_dir, atom)
        command += "> {}/POTCAR".format(path)
        sp.call(command, shell=True)
    f.close()
    print("POTCAR written")

def scf_incar(write_path, lwave=True, algo=None, prec=None, npar=None, ncore=None, 
             kpar=None, lcharg=None, lreal=None, ismear=None, sigma=None):
    '''Must be called after scf is created.'''
    if os.path.exists(f"{write_path}/POTCAR"):
        params = {
            "ALGO": algo if algo is not None else "Normal",
            "PREC": prec if prec is not None else "Accurate",
            "NPAR": npar if npar is not None else 8,
            "NCORE": ncore if ncore is not None else 8,
            "KPAR": kpar if kpar is not None else 2,
            "LREAL": lreal if lreal is not None else "Auto",
            "LWAVE": lwave if lwave is not False else True,
            "LCHARG": lcharg if lcharg is not None else False,       
            "IMSEAR": ismear if ismear is not None else 0,
            "SIGMA": sigma if sigma is not None else 0.03
        }
        
        params["NELM"] = 250
        params["NELMIN"] = 4
        params["ISTART"] = 0    
        params["EDIFF"] = 1e-7
        params["EDIFFG"] = 1e-5
        params["ISIF"] = 2 # CHANGE ME
        params["NSW"] = 0
        params["IBRION"] = -1
        params["ISYM"] = 2

        command = f"grep \"ENMAX\" {write_path}/POTCAR"
        os.system(command)
        output = os.popen(command).read();
        energy_vals = [float(i) for i in output.replace(";", "").split() if re.match(r'^-?\d+(?:\.\d+)$', i) is not None]
        params["ENCUT"] = 1.3 * max(energy_vals)

        incar = Incar.from_dict(params)
        incar.write_file(f"{write_path}/INCAR")

        print("Self-consistent field (SCF) INCAR written")

    else:
        print(f"{write_path}/POTCAR does not exist.")

def relax_incar(write_path, isif, nsw, ibrion, isym, ivdw=False, lwave=False,
              algo=None, prec=None, npar=None, ncore=None, 
             kpar=None, lreal=None, lcharg=None, ismear=None, sigma=None):
    '''Must be called after scf is created.'''
    if os.path.exists(f"{write_path}/POTCAR"):
        params = {
            "ALGO": algo if algo is not None else "Normal",
            "PREC": prec if prec is not None else "Accurate",
            "NPAR": npar if npar is not None else 8,
            "NCORE": ncore if ncore is not None else 8,
            "KPAR": kpar if kpar is not None else 2,
            "LREAL": lreal if lreal is not None else "Auto",
            "LWAVE": lwave if lwave is not False else False,
            "LCHARG": lcharg if lcharg is not None else "False",
            "IMSEAR": ismear if ismear is not None else 0,
            "SIGMA": sigma if sigma is not None else 0.03,
            "IVDW": ivdw if ivdw is not False else 0
        }
        
        params["NELM"] = 250
        params["NELMIN"] = 4
        params["ISTART"] = 0    
        params["EDIFF"] = 1e-5
        params["EDIFFG"] = -0.01
        params["ISIF"] = isif # CHANGE ME
        params["NSW"] = nsw
        params["IBRION"] = ibrion
        params["ISYM"] = isym

        command = f"grep \"ENMAX\" {write_path}/POTCAR"
        os.system(command)
        output = os.popen(command).read();
        energy_vals = [float(i) for i in output.replace(";", "").split() if re.match(r'^-?\d+(?:\.\d+)$', i) is not None]
        params["ENCUT"] = 1.3 * max(energy_vals)

        incar = Incar.from_dict(params)
        incar.write_file(f"{write_path}/INCAR")

        print("Relaxation INCAR written")

    else:
        print(f"{write_path}/POTCAR does not exist.")


def run_cpu(structure, destination, hours=18, cpu=1):
    s = f'''#!/bin/bash
    
#SBATCH -A m4064
#SBATCH -N {cpu}
#SBATCH -c 128
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t {hours}:00:00
#SBATCH -J {structure}
#SBATCH -o job.out
#SBATCH -e job.err

module load vasp/6.4.1-cpu

srun vasp_std'''
    f = open(f"{destination}/run.slurm", 'w')
    f.write(s)
    f.close()
    
def make_parchg_incar(path):#, vbm, cbm, vbm_cbm_tol):
    '''Creates an INCAR for a VASP partial charge calculation,
        specifically for creating VBM and CBM PARCHG files. It 
        searches through the OUTCAR of a ground state calculation (SCF)
        and finds the appropriate kpoint (KPUSE) and band (IBAND) indices
        that are closest to the VBM and CBM calculated from BSVasprun.
        It then writes the selected indices to an INCAR.

        path - root directory <type: str>
        bs - Pymatgen BSVasprun object <type: pymatgen.io.vasp.outputs.BSVasprun>
        vbm_cbm_tol - tolerance range to find VBM and CBM <type: list of floats>'''
    
    import re

    with open(f"{path}/scf/OUTCAR", 'r') as f:
        content = f.read()

    start_keyword = "k-point     " # start kpoint
    end_keyword = ".00000\n\n"  # end kpoint

    # Use a regular expression to find all text between the start and end keywords
    matches = re.findall(f'{re.escape(start_keyword)}(.*?){re.escape(end_keyword)}', content, re.DOTALL)

    vbm_list, cbm_list = [],[]
    ibands, kpoints = [], []

    for m in matches: # all groups of k indices and their respective band indices
        data = m.split('\n')
        for band in range(2, len(data)-1):
            if float(data[band+1].split()[2]) == 0.0 and float(data[band].split()[2]) != 0.0: #checking for band gap
                kpoints.append(int(data[0][0]))
                vbm_list.append(float(data[band].split()[1]))
                cbm_list.append(float(data[band+1].split()[1]))

                if int(data[band].split()[0]) not in ibands:
                    ibands.append(int(data[band].split()[0])) #vbm
                if int(data[band+1].split()[0]) not in ibands:
                    ibands.append(int(data[band+1].split()[0])) #vbm

                print(f"{int(data[0][0])} {data[band].split()[1]} {data[band+1].split()[1]}")

    encut = float(sp.run(f"grep \"ENCUT\" {path}/scf/OUTCAR", shell=True, capture_output=True, text=True).stdout.split()[2])
    ispin = int(sp.run(f"grep \"ISPIN\" {path}/scf/OUTCAR", shell=True, capture_output=True, text=True).stdout.split()[2])
    isym = int(sp.run(f"grep \"ISYM\" {path}/scf/OUTCAR", shell=True, capture_output=True, text=True).stdout.split()[2])

    vbm_band, cbm_band = max(vbm_list), min(cbm_list)
    vbm_index, cbm_index = vbm_list.index(vbm_band), cbm_list.index(cbm_band)
    band_gap = cbm_band - vbm_band

    print(f"VBM: {vbm_band} eV || CBM: {cbm_band} eV || band gap: {band_gap:.3f} eV")
    print(f"IBANDS: {ibands}")
    print(f"KPOINTS: {kpoints[vbm_index]} {kpoints[cbm_index]}")

    if not os.path.exists(f"{path}/parchg"):
        os.mkdir(f"{path}/parchg")
    with open(f"{path}/parchg/INCAR", 'w') as f:
        s = f"# VBM: {ibands[0]} {kpoints[vbm_index]} || CBM: {ibands[1]} {kpoints[cbm_index]}"
        s += "\nSYSTEM = STM simulation\n# must be the same as the SCF calc\n"
        s += f"\nENCUT = {encut}\nISPIN = {ispin}\nISYM = {isym}" 
        s += f"\nKPUSE = {kpoints[vbm_index]} {kpoints[cbm_index]}"
        s += f"\nIBAND = {ibands[0]} {ibands[1]}"
        s += "\n\nLPARD = .TRUE.\nLSEPB = .TRUE.\nLSEPK = .TRUE.\n"
        f.write(s)
    f.close()

    os.system(f"cp {path}/scf/POSCAR {path}/parchg/")
    print(f"PARCHG INCAR file written!")