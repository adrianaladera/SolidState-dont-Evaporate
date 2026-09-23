root = "/data/NFS/potato/aladera/COGITO/"

hl_gaps = {
    "2,6-dimethyl":   4.115,
    "1naphthyl":      2.949,
    "3methoxy":       3.861,
    "2MMB":           3.081,
    "2methoxy":       3.773,
    "2butane":        4.745,
    "2propane":       4.653,
    "gal-hydrated": 5.189,   
    # "gal-dehydrated": 4.472,
    "glu-dehydrated": 5.207,
    "glu-hydrated":   4.945,
}

from COGITO_dft.COGITO import run_cogito
from COGITO_dft.COGITOpost import run_cogito_model
import os

skipped = []
for key in hl_gaps:
        path = f"{root}/{key}/"
        if not os.path.exists(f"{path}/tb_input.txt"):
            print(f"{key} not run!")
            try:
                run_cogito(directory=f"{path}/")
            except Exception as e:
                print(f"[skip] {key}: run_cogito failed -- {e!r}")
                skipped.append(key)
                continue
        if not os.path.exists(f"{path}/all_unique_bonds.json"):
                    print(f"{key} not run!")
                    try:
                        run_cogito_model(dir=f"{path}/")
                    except Exception as e:
                        print(f"[skip] {key}: run_cogito failed -- {e!r}")
                        skipped.append(key)
                        continue

if skipped:
    print(f"\nskipped {len(skipped)} structure(s): {skipped}")