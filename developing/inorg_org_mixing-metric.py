"""
By Adriana J. Ladera
Inorganic-vs-organic PDOS mixing metric for hybrid (organic/inorganic) 1D
semiconductors, computed at the valence and conduction band edges.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun

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
          "gal-dehydrated":4.472,
          "glu-dehydrated":5.207, 
          "glu-hydrated":4.945}
COLORS = {"2,6-dimethyl":"#FF0000", 
          "1naphthyl":"#FFBE00", 
          "3methoxy":"#FFE200", 
          "2MMB":"#FCFF00", 
          "2methoxy":"#B3FF00",
          "2butane":"#77FF00",
          "2propane":"#6CFF00",
          "gal-hydrated":"#00FF38",
          "gal-dehydrated":"#00FF38",
          "glu-dehydrated":"#00FF38",            
          "glu-hydrated":"#00E2FF"}
MARKERS = {"2,6-dimethyl":"H", 
          "1naphthyl":"o", 
          "3methoxy":"v", 
          "2MMB":"*", 
          "2methoxy":"s",
          "2butane":"p",
          "2propane":"h",
          "gal-hydrated":"^",
          "gal-dehydrated":"d",
          "glu-dehydrated":"D",            
          "glu-hydrated":"8"}

INORGANIC = {"Ag", "S"}  # the 1D AgS wire (element-based default)


# #
# density summation
# #
def _sum_spins(dos):
    """Sum Spin.up (+ Spin.down if present) into a single density array."""
    return sum(dos.densities.values())


def _summed_density_by_element(edos, elements):
    """Sum spin channels and all requested elements into one density array."""
    total = None
    for el, d in edos.items():
        if el.symbol not in elements:
            continue
        dens = _sum_spins(d)
        total = dens if total is None else total + dens
    if total is None:
        raise ValueError(f"No element-projected DOS found for {elements}")
    return total


def _summed_density_by_site(complete_dos, site_indices):
    """Sum spin channels over an explicit set of site indices."""
    structure = complete_dos.structure
    total = None
    for i in site_indices:
        d = complete_dos.get_site_dos(structure[i])
        dens = _sum_spins(d)
        total = dens if total is None else total + dens
    if total is None:
        raise ValueError(f"No site-projected DOS found for sites {site_indices}")
    return total


# #
# integration
# #
def _integrate_window(energies, density, e_lo, e_hi):
    """Trapezoidal integral over [e_lo, e_hi] with interpolated endpoints, so the
    window width is identical across structures regardless of the energy grid.

    Raises if the window runs off the grid (np.interp would otherwise clamp the
    endpoint to the boundary value and silently count a region that has no data).
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


# #
# normalization
# #
def _zval(potcar_single):
    """Robust ZVAL (valence electrons per species) accessor."""
    if hasattr(potcar_single, "zval"):
        return float(potcar_single.zval)
    if hasattr(potcar_single, "nelectrons"):
        return float(potcar_single.nelectrons)
    return float(potcar_single.keywords["ZVAL"])


def _norm_factors(mode, inorganic, organic, structure, potcar):
    """Normalization denominators for the inorganic / organic partitions.

    Only supported for the element-based partition; site-based partitions should
    use normalize='atoms' (count of sites) handled by the caller, or None.
    """
    if mode is None:
        return 1.0, 1.0
    counts = structure.composition.get_el_amt_dict()
    if mode == "atoms":
        n_in = sum(v for k, v in counts.items() if k in inorganic)
        n_or = sum(v for k, v in counts.items() if k in organic)
        return n_in, n_or
    if mode == "valence":
        z = {p.element: _zval(p) for p in potcar}
        n_in = sum(counts[k] * z[k] for k in counts if k in inorganic)
        n_or = sum(counts[k] * z[k] for k in counts if k in organic)
        return n_in, n_or
    raise ValueError(f"unknown normalize mode {mode!r}")


# #
# main whatever
# #
def mixing_metric(
    complete_dos,
    window=1.0,
    vbm=None,
    cbm=None,
    inorganic=INORGANIC,
    inorganic_sites=None,   # explicit site indices -> overrides element partition
    organic_sites=None,     # optional; defaults to "all sites not inorganic"
    normalize=None,         # None (recommended) | "atoms" | "valence"
    structure=None,         # required for normalize != None (element mode)
    potcar=None,            # required for normalize == "valence"
):
    """Inorganic-vs-organic PDOS mixing metric at the valence and conduction edges.

    Partitioning

    * Default: by element symbol (`inorganic` set; organic = everything else).
      Fast, but cannot separate a shared element that lives on both sublattices.
    * Site-based: pass `inorganic_sites` (indices into complete_dos.structure).
      Use this if, e.g., sulfur appears in both the wire and an organic linker.

    Returns
-
    {"valence": {...}, "conduction": {...}} with, per edge:
        E_window            (e_lo, e_hi) in absolute eV
        I_inorganic         integrated inorganic PDOS
        I_organic           integrated organic PDOS
        ratio_in_over_or    I_in / I_or  (0 => organic band, inf => inorganic band)
        frac_inorganic      I_in / (I_in + I_or)  in [0, 1]
        mixing_balance      2*sqrt(I_in*I_or)/(I_in+I_or) in [0, 1]
                            (0 => pure/unmixed band, 1 => equal contribution)
    """
    energies = np.asarray(complete_dos.energies)  # absolute eV, matches get_cbm_vbm

    if vbm is None or cbm is None:
        cbm_auto, vbm_auto = complete_dos.get_cbm_vbm()  # returns (cbm, vbm)
        vbm = vbm_auto if vbm is None else vbm
        cbm = cbm_auto if cbm is None else cbm

    # build the two summed densities--
    if inorganic_sites is not None:
        inorganic_sites = list(inorganic_sites)
        n_sites = len(complete_dos.structure)
        if organic_sites is None:
            organic_sites = [i for i in range(n_sites) if i not in set(inorganic_sites)]
        else:
            organic_sites = list(organic_sites)
        d_in = _summed_density_by_site(complete_dos, inorganic_sites)
        d_or = _summed_density_by_site(complete_dos, organic_sites)
        if normalize == "atoms":
            d_in = d_in / len(inorganic_sites)
            d_or = d_or / len(organic_sites)
        elif normalize == "valence":
            raise ValueError(
                "normalize='valence' is only implemented for element-based "
                "partitioning; use normalize='atoms' or None with sites."
            )
    else:
        edos = complete_dos.get_element_dos()            # {Element: Dos}
        organic = {el.symbol for el in edos} - set(inorganic)
        d_in = _summed_density_by_element(edos, inorganic)
        d_or = _summed_density_by_element(edos, organic)
        n_in, n_or = _norm_factors(normalize, inorganic, organic, structure, potcar)
        d_in, d_or = d_in / n_in, d_or / n_or

    # integrate over each edge window-
    windows = {
        "valence":    (vbm - window, vbm), # within 1eV of VBM
       "conduction": (cbm, cbm + window), # within 1eV of CBM
    }
    out = {}
    for name, (e_lo, e_hi) in windows.items():
        I_in = _integrate_window(energies, d_in, e_lo, e_hi)
        I_or = _integrate_window(energies, d_or, e_lo, e_hi)
        tot = I_in + I_or
        out[name] = {
            "E_window": (float(e_lo), float(e_hi)),
            "I_inorganic": I_in,
            "I_organic": I_or,
            "ratio_in_over_or": (I_in / I_or) if I_or else np.inf,
            "frac_inorganic": (I_in / tot) if tot else np.nan,
            "mixing_balance": (2.0 * np.sqrt(I_in * I_or) / tot) if tot else np.nan,
        }
    return out


if __name__ == "__main__":

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    edge_axes = {"valence": axes[0], "conduction": axes[1]}
    edge_data = {"valence": {"x": [], "y": []}, "conduction": {"x": [], "y": []}}

    for key, gap in hl_gaps.items():
        vasprun_path = os.path.join(root, key, "dos", "vasprun.xml")
        if not os.path.exists(vasprun_path):
            print(f"skipping {key!r}: no vasprun.xml found at {vasprun_path}")
            continue

        vr = Vasprun(vasprun_path)
        cdos = vr.complete_dos
        result = mixing_metric(cdos)
        # print(f"{key}: {result}")

        for edge, ax in edge_axes.items():
            ratio = result[edge]["ratio_in_over_or"]
            edge_data[edge]["x"].append(gap)
            edge_data[edge]["y"].append(ratio)
            ax.scatter(
                gap,
                ratio,
                color=COLORS[key],
                marker=MARKERS[key],
                s=80,
                edgecolors="black",
                linewidths=0.5,
                label=key,
            )

    for edge, ax in edge_axes.items():
        ax.set_title(f"{edge.capitalize()} band edge")
        ax.set_xlabel("HOMO-LUMO gap (eV)")
        ax.set_ylabel("I$_{inorganic}$ / I$_{organic}$")

        x, y = edge_data[edge]["x"], edge_data[edge]["y"]
        if len(x) >= 2:
            r = np.corrcoef(x, y)[0, 1]
            ax.text(
                0.95, 0.95,
                f"$r$ = {r:.3f}\n$r^2$ = {r**2:.3f}",
                transform=ax.transAxes,
                ha="right", va="top",
            )

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), title="Structure")
    fig.tight_layout()
    fig.savefig(f"inorg_org_mixing_metric{TAG}.png", dpi=300, bbox_inches="tight")
    plt.show()