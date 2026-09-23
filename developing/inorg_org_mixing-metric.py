"""
By Adriana J. Ladera
Inorganic-vs-organic PDOS mixing metric for hybrid (organic/inorganic) 1D
semiconductors, computed at the valence and conduction band edges.
"""

import functools
import numpy as np
import os
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.dos import add_densities
from COGITO_dft.COGITO import run_cogito
from COGITO_dft.COGITOpost import run_cogito_model
from COGITO_dft.COGITOpost import COGITO_TB_Model as CoTB, COGITO_UNIFORM

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
TAG = "_no-outlier"

root = "/data/NFS/potato/aladera/1D_rainbow/"
hl_gaps = {"2,6-dimethyl":4.115, 
          "1naphthyl":2.949, 
          "3methoxy":3.861, 
          "2MMB":3.081, 
          "2methoxy":3.773,
          "2butane":4.745,
          "2propane":4.653,
          "gal-hydrated":5.189, # OUTLIER
        #   "gal-dehydrated":4.472,
        #   "glu-dehydrated":5.207, 
          "glu-hydrated":4.945}
COLORS = {"2,6-dimethyl":"#FF0000", 
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
MARKERS = {"2,6-dimethyl":"H", 
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

INORGANIC = {"Ag", "S"}  # the 1D AgS wire (element-based default)


# density summation
def _sum_spins(dos):
    """Sum Spin.up (+ Spin.down if present) into a single density array."""
    return sum(dos.densities.values())

# integration
def _integrate_window(energies, density, e_lo, e_hi):
    """Trapezoidal integral over [e_lo, e_hi] with interpolated endpoints, so the
    window width is identical across structures regardless of the energy grid.
    """
    if e_lo < energies[0] or e_hi > energies[-1]:
        raise ValueError(
            f"window [{e_lo:.3f}, {e_hi:.3f}] eV exceeds DOS energy range "
            f"[{energies[0]:.3f}, {energies[-1]:.3f}] eV; increase the energy "
            f"range in your DOS calc or shrink the window."
        )
    # window bounds for either VBM or CBM
    d_lo = np.interp(e_lo, energies, density)
    d_hi = np.interp(e_hi, energies, density)
    mask = (energies > e_lo) & (energies < e_hi)
    e = np.concatenate(([e_lo], energies[mask], [e_hi]))
    d = np.concatenate(([d_lo], density[mask], [d_hi]))
    return float(_trapz(d, e))


# normalization
def _zval(potcar_single):
    """Robust ZVAL (valence electrons per species) accessor."""
    if hasattr(potcar_single, "zval"):
        return float(potcar_single.zval)
    if hasattr(potcar_single, "nelectrons"):
        return float(potcar_single.nelectrons)
    return float(potcar_single.keywords["ZVAL"])


# smooth band-edge weighting (replaces the hard 1 eV step window)
def _edge_weight(energies, edge, kind, width=1.0, shape="tanh"):
    """Smooth 'turn-off' weight in [0,1] for states near a band edge.
        Per Emoomoo's rec: a step function makes the result jump when the
    cutoff crosses a band that starts ~1 eV in (like 2,6-dimethyl VB). 
    Smooth turn-off also physically represents deeper/higher
    states being progressively less accessible.

    energies : absolute-eV grid (same frame as `edge`)
    edge     : VBM or CBM in absolute eV
    kind     : "valence" -> weight states at/below VBM, decaying downward
               "conduction" -> weight states at/above CBM, decaying upward
    width    : decay length in eV (~ your old 1.0 eV hard cutoff)
    shape    : "tanh" (soft shoulder ~1 to ~width then rolls to 0) or
               "exp"  (1 at the edge, exp(-d/width) into the band)
    """
    e = np.asarray(energies, dtype=float)
    if kind == "valence":
        d = edge - e # >0 below the VBM
    elif kind == "conduction":
        d = e - edge # >0 above the CBM
    else:
        raise ValueError("kind must be 'valence' or 'conduction'")

    if shape == "exp":
        w = np.where(d >= 0.0, np.exp(-d / width), 0.0)
    elif shape == "tanh":
        w = 0.5 * (1.0 - np.tanh((d - width) / (0.5 * width)))
        w = np.where(d >= -0.25 * width, w, 0.0)   # kill wrong-side leakage
    else:
        raise ValueError("Sorry brother bear, shape must be 'tanh' or 'exp'")
    return w


def _weighted_integral(energies, density, edge, kind, width=1.0, shape="tanh"):
    """Trapezoidal integral of density times the smooth edge weight."""
    e = np.asarray(energies, dtype=float)
    order = np.argsort(e)
    e_s = e[order]
    y_s = (np.asarray(density, dtype=float) * _edge_weight(e, edge, kind, width, shape))[order]
    return float(_trapz(y_s, e_s))


# density summation
def _sum_spins(dos):
    """Sum Spin.up (+ Spin.down if present) into one density array."""
    return sum(dos.densities.values())


def partition_sites(structure, inorganic_elements=INORGANIC,
                    interface_element="S", interface_to="inorganic"):
    """Split site indices into inorganic vs organic to explicitly list 
        where you want the S to be counted.

    interface_to : "organic" | "inorganic" | "split"
        "organic"   -> the C-bonded S is treated as part of the ligand fragment
        "inorganic" -> S stays with the AgS wire (your old element-default behavior)
        "split"     -> S excluded from both (report only Ag vs C/N/O/H)

    Returns (inorganic_sites, organic_sites) as lists of indices.
    """
    inorg, org = [], []
    for i, site in enumerate(structure):
        sym = site.specie.symbol if hasattr(site, "specie") else site.species_string
        if sym == interface_element:
            if interface_to == "inorganic":
                inorg.append(i)
            elif interface_to == "organic":
                org.append(i)
            elif interface_to == "split":
                continue
            else:
                raise ValueError(f"bad interface_to={interface_to!r}")
        elif sym in inorganic_elements:
            inorg.append(i)
        else:
            org.append(i)
    return inorg, org


def _summed_density_by_site(complete_dos, site_indices):
    """Sum spin channels over an explicit set of site indices.
    """
    structure = complete_dos.structure
    total = None
    for i in site_indices:
        site = structure[i]
        if site in complete_dos.pdos:
            dens = _sum_spins(complete_dos.get_site_dos(site))
        else:
            site_pdos = complete_dos.pdos[i]
            per_orbital = functools.reduce(add_densities, site_pdos.values())
            dens = _sum_spins_dict(per_orbital)
        total = dens if total is None else total + dens
    if total is None:
        raise ValueError(f"no site-projected DOS for sites {site_indices}")
    return total


def _sum_spins_dict(densities):
    """Sum Spin.up (+ Spin.down if present) from a raw {Spin: array} dict."""
    return sum(densities.values())


# main whatever
def mixing_metric(
    complete_dos,
    window=1.0,
    shape="tanh", # "tanh" | "exp" | "step" (the OG way w/hard 1eV cutoff)
    vbm=None,
    cbm=None,
    inorganic_elements=INORGANIC,
    interface_element="S",
    interface_to="organic",    # where the thiolate S goes; see partition_sites
    inorganic_sites=None,      # fill in if you wanna override the interface bs
    organic_sites=None,
    normalize="atoms",         # "atoms" | None  (cancels in the ratio; matters for %)
):
    """Inorganic-vs-organic DOS mixing at the valence and conduction edges.


        E_edge            VBM or CBM (absolute eV)
        I_inorganic       smooth-windowed integral of the inorganic DOS
        I_organic         smooth-windowed integral of the organic DOS
        pct_organic       100 * I_or / (I_or + I_in)      -> [0, 100]
        pct_inorganic     100 - pct_organic               -> [0, 100]
        ratio_in_over_or  I_in / I_or  (legacy metric; unchanged meaning)
    """
    energies = np.asarray(complete_dos.energies)   # absolute eV
    structure = complete_dos.structure

    if vbm is None or cbm is None:
        cbm_auto, vbm_auto = complete_dos.get_cbm_vbm()   # (cbm, vbm)
        vbm = vbm_auto if vbm is None else vbm
        cbm = cbm_auto if cbm is None else cbm

    # choose site partition
    if inorganic_sites is None:
        inorganic_sites, organic_sites = partition_sites(
            structure, inorganic_elements, interface_element, interface_to
        )
    elif organic_sites is None:
        allset = set(range(len(structure)))
        organic_sites = sorted(allset - set(inorganic_sites))

    d_in = _summed_density_by_site(complete_dos, inorganic_sites)
    d_or = _summed_density_by_site(complete_dos, organic_sites)

    if normalize == "atoms":
        # cancels in the ratio, but makes I_in / I_or comparable in the % metric
        d_in = d_in / max(len(inorganic_sites), 1)
        d_or = d_or / max(len(organic_sites), 1)
    elif normalize is not None:
        raise ValueError("normalize must be 'atoms' or None with a site partition")

    edges = {"valence": ("valence", vbm), "conduction": ("conduction", cbm)}
    out = {}
    for name, (kind, edge) in edges.items():
        if shape == "step":   # OG version with hard 1eV
            e_lo, e_hi = (edge - window, edge) if kind == "valence" else (edge, edge + window)
            I_in = _integrate_window(energies, d_in, e_lo, e_hi)
            I_or = _integrate_window(energies, d_or, e_lo, e_hi)
        else:
            I_in = _weighted_integral(energies, d_in, edge, kind, window, shape)
            I_or = _weighted_integral(energies, d_or, edge, kind, window, shape)
        tot = I_in + I_or
        out[name] = {
            "E_edge": float(edge),
            "I_inorganic": I_in,
            "I_organic": I_or,
            "pct_organic": (100.0 * I_or / tot) if tot else np.nan,
            "pct_inorganic": (100.0 * I_in / tot) if tot else np.nan,
            "ratio_in_over_or": (I_in / I_or) if I_or else np.inf,
        }
    return out

if __name__ == "__main__":
    from pymatgen.io.vasp.outputs import Vasprun  # keep for the VASP path

    USE_COGITO = True     # flip to compare VASP vs COGITO projections
    METRIC_KEY = "pct_organic"   # or "ratio_in_over_or" for the legacy plot

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    edge_axes = {"valence": axes[0], "conduction": axes[1]}
    edge_data = {"valence": {"x": [], "y": []}, "conduction": {"x": [], "y": []}}

    for key, gap in hl_gaps.items():
        if USE_COGITO:
            path = f"{root}/{key}"
            if not os.path.exists(f"{path}/tb_input.txt"):
                print(f"{key} not run!")
                run_cogito(directory=path)
                run_cogito_model(dir=path)
            
            my_CoTB = CoTB(path)
            my_CoTB.normalize_params()                      # precise normalization
            my_CoTB.restrict_params(maximum_dist=15, minimum_value=0.00001)
            
                    # 2. choose a k-grid: denser along the 1D chain axis
                    # k-density should scale like 1/|a_i|; the chain axis (smallest |a|) gets
                    # the most points. This prints a suggestion — set GRID explicitly and then
                    # densify until the edge-window integrals stop moving.
            lengths = np.linalg.norm(my_CoTB._a, axis=1) # getting along shortest lat vec || to 1D chain
            base = 24.0                                      # tune: higher = denser
            suggested = tuple(max(2, int(round(base / L))) for L in lengths)
            print(f"lattice |a_i| = {np.round(lengths, 3)}  ->  suggested GRID = {suggested}")
            
            # build realspace grid
            GRID = suggested # or hard-code to some bullshit like (2, 2, 10)
            uni = COGITO_UNIFORM(my_CoTB, GRID)
            cdos = uni.get_pymatgen_completedos()

        else:
            vasprun_path = os.path.join(root, key, "dos", "vasprun.xml")
            if not os.path.exists(vasprun_path):
                print(f"skipping {key!r}: no vasprun.xml at {vasprun_path}")
                continue
            cdos = Vasprun(vasprun_path).complete_dos

        result = mixing_metric(cdos, window=1.0, shape="tanh",
                               interface_to="organic", normalize="atoms")

        for edge, ax in edge_axes.items():
            val = result[edge][METRIC_KEY]
            edge_data[edge]["x"].append(gap)
            edge_data[edge]["y"].append(val)
            ax.scatter(gap, val, color=COLORS[key], marker=MARKERS[key],
                       s=80, edgecolors="black", linewidths=0.5, label=key)

    ylabel = ("% organic character at edge" if METRIC_KEY == "pct_organic"
              else "I$_{inorganic}$ / I$_{organic}$")
    for edge, ax in edge_axes.items():
        ax.set_title(f"{edge.capitalize()} band edge")
        ax.set_xlabel("HOMO-LUMO gap (eV)")
        ax.set_ylabel(ylabel)
        x, y = edge_data[edge]["x"], edge_data[edge]["y"]
        if len(x) >= 2:
            r = np.corrcoef(x, y)[0, 1]
            ax.text(0.95, 0.95, f"$r$ = {r:.3f}\n$r^2$ = {r**2:.3f}",
                    transform=ax.transAxes, ha="right", va="top")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left",
               bbox_to_anchor=(1.0, 0.5), title="Structure")
    fig.tight_layout()
    fig.savefig(f"inorg_org_mixing_metric{TAG}.png", dpi=300, bbox_inches="tight")
    plt.show()