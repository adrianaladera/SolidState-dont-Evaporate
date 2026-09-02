from pymatgen.io.vasp.outputs import BSVasprun
import joblib
import hashlib
import os

METALS = ['Ag', 'Au', 'Hg', 'Cu']
CHALCS = ['Te', 'Se', "S"]

def sort_elements(poscar):
    '''Sorts elements between inorganic elements (metals, chalcogens)
        and organic elements.
        poscar: POSCAR file with VASP input'''
    inorganics = METALS + CHALCS
    inorgs, orgs = [], []
    with open(poscar) as f:
        lines = f.readlines()
        cunt = 0
        flag = 0
        for line in lines:
            for i in line.split():
                if str(i) in inorganics and cunt > 1:
                    flag = 1
            if flag:
                break
            cunt += 1
        for i in lines[cunt].split():
            if i in inorganics:
                inorgs.append(i)
            else:
                orgs.append(i)
    return inorgs, orgs

def _get_cache_path(xml_path: str, cache_dir: str = ".cache") -> str:
    mtime = os.path.getmtime(xml_path)
    key = hashlib.md5(f"{xml_path}:{mtime}".encode()).hexdigest()
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, f"{key}.pkl")

def load_band_structure_cached(xml_path: str, kpoints_path: str, **kwargs):
    cache_path = _get_cache_path(xml_path)
    if os.path.exists(cache_path):
        return joblib.load(cache_path)
    
    bandrun = BSVasprun(xml_path, **kwargs)
    bs = bandrun.get_band_structure(kpoints_path)
    joblib.dump(bs, cache_path)
    return bs