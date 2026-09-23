#!/usr/bin/env python3
"""
Extract ONE scalar C-S bond-strength descriptor per MOCha.

For each structure we read COGITO's `all_unique_bonds.json` (made by
run_cogito_model / `COGITOpost.py --dir DIR`). Every entry is a
symmetry-unique bond:

    key   = "<elem1> <elem2> <dist> <deg>"      e.g. "C S 1.81 2.0"
    value = { "cohp": <iCOHP, eV per bond>,
              "coop": <iCOOP, elec per bond>,
              "degeneracy": <# of this bond per cell>,
              "bond length": <Angstrom>,
              "all bonds": [...] }

iCOHP sign convention (COGITO / standard COHP):
    negative = bonding,  positive = antibonding
=> a MORE NEGATIVE number means a STRONGER C-S bond.

We reduce the (possibly several) C-S entries down to a single number three ways;
pick whichever is the right descriptor for ligma correlation:

    reduce="sum"  -> total C-S bonding in the cell   = sum_i cohp_i * deg_i
    reduce="mean" -> degeneracy-weighted mean per bond = sum_i cohp_i*deg_i / sum_i deg_i
    reduce="nn"   -> the single nearest-neighbour (shortest) C-S bond, per bond

Set SIGNED=False to flip the sign so that "bigger = stronger" for plotting.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt


# data location + HOMO-LUMO gaps
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
    # "glu-dehydrated": 5.207,
    "glu-hydrated":   4.945,
}
COLORS = {
        "2,6-dimethyl":"#FF0000", 
          "1naphthyl":"#FFBE00", 
          "3methoxy":"#FFE200", 
          "2MMB":"#FCFF00", 
          "2methoxy":"#B3FF00",
          "2butane":"#77FF00",
          "2propane":"#6CFF00",
          "gal-hydrated":"#00FF38",
        #   "gal-dehydrated":"#00FF38",
        #   "glu-dehydrated":"#00FF38",            
          "glu-hydrated":"#00E2FF"}
MARKERS = {
        "2,6-dimethyl":"H", 
          "1naphthyl":"o", 
          "3methoxy":"v", 
          "2MMB":"*", 
          "2methoxy":"s",
          "2butane":"p",
          "2propane":"h",
          "gal-hydrated":"^",
        #   "gal-dehydrated":"d",
        #   "glu-dehydrated":"D",            
          "glu-hydrated":"8"}

# tuning shit
REDUCE = "nn"    # "sum" | "mean" | "nn"
SIGNED = True     # True -> negative = bonding ;  False -> positive = stronger
GENERATE_IF_MISSING = False   # run COGITO to build the json if it isn't there yet


# core: read one all_unique_bonds.json -> scalar C-S bond strength
def cs_bond_strength(bond_json, a="C", b="S", reduce="nn", signed=True,
                     verbose=False):
    """Return a single scalar describing the a-b (default C-S) bond strength 
     or np.nan if the structure has no a-b bond in the file.
    
        a: atom A
        b: atom B
        reduce: "nn" nearest neigbussy, "sum" total a-b bonding per unit cell,
                "mean" degeneracy-weighted mean per bond
    """
    with open(bond_json, "r") as f:
        bonds = json.load(f)

    dists, cohps, degs = [], [], []
    for name, info in bonds.items():
        parts = name.split()
        e1, e2 = parts[0], parts[1]
        if {e1, e2} != {a, b}:
            continue

        cohp = info["cohp"]                       # iCOHP, eV per bond
        deg = info.get("degeneracy", 1.0)
        dist = info.get("bond length",
                        float(parts[2]) if len(parts) > 2 else np.nan)
        dists.append(dist)
        cohps.append(cohp)
        degs.append(deg)

    if not cohps:
        return np.nan

    dists = np.asarray(dists, dtype=float)
    cohps = np.asarray(cohps, dtype=float)
    degs = np.asarray(degs, dtype=float)

    if verbose:
        for d, c, g in sorted(zip(dists, cohps, degs)):
            print(f"      {a}-{b}  d={d:6.2f} A   iCOHP={c:+.4f} eV/bond   x{g:g}")

    if reduce == "sum":       # total a-b bonding per unit cell
        val = np.sum(cohps * degs)
    elif reduce == "mean":    # degeneracy-weighted mean per bond
        val = np.sum(cohps * degs) / np.sum(degs)
    elif reduce == "nn":      # nearest-neighbour (shortest) bond only, per bond
        finite = np.isfinite(dists)
        if not finite.any():
            return np.nan
        val = cohps[finite][np.argmin(dists[finite])]
    else:
        raise ValueError(f"unknown reduce={reduce!r}")

    return val if signed else -val


# optional func to build all_unique_bonds.json via COGITO if it's missing
def ensure_bond_json(dir_, **cogito_kwargs):
    path = os.path.join(dir_, "all_unique_bonds.json")
    if os.path.isfile(path):
        return path
    from COGITO_dft.COGITOpost import run_cogito_model
    d = dir_ if dir_.endswith("/") else dir_ + "/"
    run_cogito_model(dir=d, save_crystal_bonds=True, **cogito_kwargs)
    return path


if __name__ == "__main__":
    labels, gaps, cs = [], [], []

    for key, gap in hl_gaps.items():
        d = os.path.join(root, key)
        path = os.path.join(d, "all_unique_bonds.json")

        if not os.path.isfile(path):
            if GENERATE_IF_MISSING:
                path = ensure_bond_json(d)
            else:
                print(f"[skip] {key}: {path} not found")
                continue

        print(f"{key}:")
        val = cs_bond_strength(path, a="C", b="S",
                               reduce=REDUCE, signed=SIGNED, verbose=True)

        labels.append(key)
        gaps.append(gap)
        cs.append(val)
        print(f"   -> C-S ({REDUCE}) = {val:+.4f} eV\n")

    gaps = np.asarray(gaps, dtype=float)
    cs = np.asarray(cs, dtype=float)

    # dict:  {ligma: scalar}
    cs_by_ligand = dict(zip(labels, cs))
    print("cs_by_ligand =", {k: round(float(v), 4) for k, v in cs_by_ligand.items()})

    # quick correlation with the gap (drop any NaNs first)
    mask = ~np.isnan(cs)
    if mask.sum() >= 2:
        r = np.corrcoef(gaps[mask], cs[mask])[0, 1]
        print(f"\nPearson r (gap vs C-S {REDUCE}) = {r:+.3f}")


    # plot C-S bond strength vs HOMO-LUMO gap

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    for x, y, lab in zip(gaps, cs, labels):
        if np.isnan(y):
            continue
        ax.scatter(x, y, color=COLORS[lab], marker=MARKERS[lab], s=80,
                   edgecolors="black", linewidths=0.8, label=lab)
        ax.annotate(lab, (x, y), fontsize=8,
                    textcoords="offset points", xytext=(4, 4))

    if mask.sum() >= 2:
        ax.text(0.05, 0.95, f"$R$ = {r:.3f}\n$R^2$ = {r**2:.3f}",
                transform=ax.transAxes, ha="left", va="top", fontsize="small",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    sign_note = "neg = stronger bonding" if SIGNED else "pos = stronger bonding"
    ax.set_xlabel("HOMO-LUMO gap (eV)")
    ax.set_ylabel(f"C-S iCOHP  [{REDUCE}]  (eV)   ({sign_note})")
    ax.set_title("C-S bond strength vs HOMO-LUMO gap")
    # ax.axhline(0, linewidth=0.8, color="k")
    # ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig("CS_bond_vs_gap.png", dpi=300, bbox_inches="tight")
    plt.show()