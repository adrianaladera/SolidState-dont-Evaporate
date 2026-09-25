''' COHP_COOP.py
By Adriana J. Ladera

Takes the outputs produced from a VASP run and passes
theme through the COGITO atomic orbital basis. Obtains
the valence and conduction band edges, and then integrates
the COHP within a window of 1 eV from the VBM and 1 eV from 
the CBM. Plots the integrated COHP value against the HOMO-LUMO
gaps of the protonated thiol ligand (gas phase) to produce 
correlation plots for the valence and conduction bands, saved
as "cohp_vs_hlgap.png".

'''

import os
import numpy as np
import matplotlib.pyplot as plt
from COGITO_dft.COGITO import run_cogito
from COGITO_dft.COGITOpost import run_cogito_model
from COGITO_dft.COGITOpost import COGITO_TB_Model as CoTB, COGITO_UNIFORM

root = "/data/norm_fucktorS/potato/aladera/COGITO/"
elec_band_gaps = {"2,6-dimethyl":1.546, 
          "1naphthyl":1.735, 
          "3methoxy":1.938, 
          "2MMB":1.780, 
          "2methoxy":1.832,
          "2butane":2.175,
          "2propane":2.140,
          "gal-hydrated":2.832,
        #   "gal-dehydrated":2.453,
        #   "glu-dehydrated":2.478, 
          "glu-hydrated":2.440}
hl_gaps = {"2,6-dimethyl":4.115, 
          "1naphthyl":2.949, 
          "3methoxy":3.861, 
          "2MMB":3.081, 
          "2methoxy":3.773,
          "2butane":4.745,
          "2propane":4.653,
          "gal-hydrated":5.189,
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
_trapz = getattr(np, "trapezoid", None) or np.trapz

# fragments: AgS wire vs everything on the ligand
INORG = {"Ag": ["s", "p", "d"], "S": ["s", "p"]}
ORG   = {"C": ["s", "p"], "N": ["s", "p"], "O": ["s", "p"], "H": ["s"]}

def resolve_fermi(uni):
    """E_F in the SAME frame as uni.eigvals. COGITO stores eigvals shifted by
    energy_shift; uni.efermi is in the unshifted frame, so add the shift.
    Basically gets the same reference frame for uni as COGITO eigs frame"""
    shift = float(getattr(uni, "energy_shift", 0.0) or 0.0)
    return float(uni.efermi) + shift, shift


def band_edges_robust(uni, expected_gap=None, tol=1e-6):
    """VBM/CBM on the E_F=0 axis that get_COHP_DOS uses. Cross-checks the
    Fermi-split edges against band counting and guards the fake-tiny-gap case."""

    # need eigvals bc band edges are extremal eigvals
    eigs = np.asarray(uni.eigvals, dtype=float)          # (nspin, nkpt, nband)
    eigs2d = (np.concatenate([eigs[s] for s in range(eigs.shape[0])], axis=1) # format array bullshit (3D--> 2D)
              if eigs.ndim == 3 else eigs)
    ef, shift = resolve_fermi(uni) # get E_F in same ref frame as uni.eigvals

    n_occ_k   = np.sum(eigs2d <= ef + tol, axis=1) # per-k occupied count
    n_occ_mean = float(np.mean(n_occ_k)) # diagnostic (may be non-integer)
    n_occ      = int(round(n_occ_mean)) # integer for band-counting index
    constant   = bool(np.all(n_occ_k == n_occ))
    print(f"  n_occ = {n_occ} (mean {n_occ_mean:.3f}; constant across k: {constant})")

    # split by max occupied set for VBM and min unoccupied cunt set for CBM
    occ, unocc = eigs2d[eigs2d <= ef + tol], eigs2d[eigs2d > ef + tol]
    vbm_abs, cbm_abs = float(occ.max()), float(unocc.min())     # Fermi-split
    # or band counting to get BG w/o ref to ef (just needs occupied count)
    # s = np.sort(eigs2d, axis=1)
    # vbm_bc, cbm_bc = float(s[:, n_occ - 1].max()), float(s[:, n_occ].min())  # count

    vbm_rel, cbm_rel = vbm_abs - ef, cbm_abs - ef # shift
    gap = cbm_rel - vbm_rel

    print(f"  eigvals range : [{eigs2d.min():.3f}, {eigs2d.max():.3f}] (shifted)")
    print(f"  E_F (eig frame) = {ef:.4f}  [uni.efermi={uni.efermi:.4f} + shift={shift:.4f}]")
    print(f"  n_occ = {n_occ}  (constant across k: {constant})")
    # print(f"  gap: fermi-split = {gap:.4f} eV | band-count = {cbm_bc - vbm_bc:.4f} eV")

    # check to see if this shit matches the actual electronic band gap from DFT
    if not constant:
        print("  WARNING: occupied count varies across k — metallic or E_F still off.")
    assert gap > 0.05, (f"Gap {gap:.4f} eV ~ 0 -> E_F in wrong frame; check energy_shift.")
    if expected_gap is not None and abs(gap - expected_gap) > 0.1:
        print(f"  WARNING: gap {gap:.4f} eV != expected {expected_gap} eV.")
    return vbm_rel, cbm_rel

def _weighted_integral(energies, y, e_lo, e_hi, edge, kind, width=1.0, shape="tanh"):
    """Trapezoid over [e_lo, e_hi] with interpolated endpoints so the
    window is exactly `window` eV wide on every structure's grid.
    Sort to get correct energy order of np.interp and np.trapz, otherwise
    integral is lowk useless. Interpolate since window endpoints don't 
    always exactly land on grid points, which would vary slightly per struct
    and bias the correlation metric!!"""
    order = np.argsort(energies) #
    E, Y = np.asarray(energies)[order], np.asarray(y)[order] # maps energy-y pairing even after sorting 
    y_lo, y_hi = np.interp(e_lo, E, Y), np.interp(e_hi, E, Y)
    m = (E > e_lo) & (E < e_hi) # mask to avoid double-counting grid point
    xs = np.concatenate(([e_lo], E[m], [e_hi]))
    ys = np.concatenate(([y_lo], Y[m], [y_hi])) 

    w = _edge_weight(xs, edge, kind, width, shape)
    w_integral = float(_trapz(ys*w, xs))
    w_average = float(_trapz(w, xs)) # idk exactly how to count the # of points in a DOS; in COGITO
    # it returns the DOS per energy level but technically there can be multiple states in a single energy
    # level no? in a bandstructure there are multiple k points and bands, but in DOS its just the density
    # of states at that energy? i only did this because I remember Emoo saying something about also applying 
    # the weight to the normalization...right now I'm kinda just interpreting normalization's "total number
    # of points within window" as individual energies corresponding to each DOS point. Prolly wrong lol.
    
    return w_integral, w_average

def _integrate_window(energies, y, e_lo, e_hi):
    """Trapezoid over [e_lo, e_hi] with interpolated endpoints so the
    window is exactly `window` eV wide on every structure's grid.
    Sort to get correct energy order of np.interp and np.trapz, otherwise
    integral is lowk useless. Interpolate since window endpoints don't 
    always exactly land on grid points, which would vary slightly per struct
    and bias the correlation metric!!"""
    order = np.argsort(energies) #
    E, Y = np.asarray(energies)[order], np.asarray(y)[order] # maps energy-y pairing even after sorting 
    y_lo, y_hi = np.interp(e_lo, E, Y), np.interp(e_hi, E, Y)
    m = (E > e_lo) & (E < e_hi) # mask to avoid double-counting grid point
    xs = np.concatenate(([e_lo], E[m], [e_hi]))
    ys = np.concatenate(([y_lo], Y[m], [y_hi])) 

    return float(_trapz(ys, xs))


def cohp_norm_factor(uni, mode="window", window_den=None,
                     interface_pair=("C", "S"), sigma=0.05, _struct_cache={}):
    """Normalize this hoe by number of S atoms. Fuck interface pair
        bc the C-S bond lengths vary per MOCha and I don't wanna risk
        undercounting, besides num C-S bonds == num S atoms.
   mode = "window" : weight-sum integrating w*E within the window  -> turns the integral
                       into a weighted mean (what 'normalize by # densities in
                       the window' means??)
           "s_atoms" : number of B-element (S) atoms
           "atoms" : total atoms  (Emily said prob not good wahh)
           None : 1.0, un-normalizied rawdogging
    """
    if mode is None:
        return 1.0

    if mode == "window":
        if window_den is None:
            raise ValueError("mode='window' needs window_den (the weight-sum ∫w dE)")
        return float(window_den) if abs(window_den) > 1e-12 else 1.0

    # site-count modes need the structure
    key = id(uni)
    struct = _struct_cache.get(key)
    if struct is None:
        cdos = COGITO_UNIFORM.get_pymatgen_completedos(uni, sigma=sigma, mulliken=True)
        struct = cdos.structure
        _struct_cache[key] = struct

    _, bob = interface_pair

    def _sym(site):
        return site.specie.symbol if hasattr(site, "specie") else site.species_string

    if mode == "atoms":
        return float(len(struct))
    if mode == "s_atoms":
        n = sum(1 for s in struct if _sym(s) == bob)
        return float(max(n, 1))

    raise ValueError(f"unknown norm mode {mode!r}")

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

    if shape == "exp": # exp decay further from the band edge
        # d=0 -> exp(0)=1, full weight, less further from band edge
        w = np.where(d >= 0.0, np.exp(-d / width), 0.0)
    elif shape == "tanh":
        # tan range [-1,1] --> 1-tanh range [2,0] * 0.5 = [1,0]
        w = 0.5 * (1.0 - np.tanh((d - width) / (0.5 * width)))
        # since tan never really 0/1, .042069 arbitrary hard cutoff for 
        # 0.042069 beyond width, then force 0
        w = np.where(d >= -0.042069 * width, w, 0.0)   # kill wrong-side leakage
    else:
        raise ValueError("Sorry brother bear, shape must be 'tanh' or 'exp'")
    return w

def cross_cohp_window(uni, vbm_rel, cbm_rel, window=1.0,
                      inorg=INORG, org=ORG,
                      max_dist=3.2, sigma=0.05, spin=0,
                      norm_mode="window", shape="tanh"):
    """Inorganic-organic COHP integrated over the valence and conduction
    edge windows. Huzz are normalized so MOChas are comparable.
    """

    # use COGITO_UNIFORM COHP to get COHP stuff
    energies, cohp_by_spin = COGITO_UNIFORM.get_COHP_DOS(
        uni, orbs=[inorg, org], max_dist=max_dist, sigma=sigma,
        include_onsite=False, save_plot=False)
    cohp = np.asarray(cohp_by_spin[spin])

    windows = {"valence": ("valence", vbm_rel - window, vbm_rel),
               "conduction": ("conduction", cbm_rel, cbm_rel + window)}
    edges = {"valence": vbm_rel, "conduction": cbm_rel}

    im_coming_out = {}
    for name, (kind, lo, hi) in windows.items():
        w_integral, w_avg = _weighted_integral(
            energies, cohp, lo, hi, edges[name], kind, window, shape)

        if norm_mode == "window":
            norm_fucktor = cohp_norm_factor(uni, mode="window", window_den=w_avg)
            # eventually migrate normalization calculation into cohp_norm_factor()
            # rather than calculating it within the integral for consistency
        else:
            norm_fucktor = cohp_norm_factor(uni, mode=norm_mode,
                                  interface_pair=("C", "S"), sigma=sigma)

        im_coming_out[name] = {"window": (float(lo), float(hi)),
                     "int_COHP": w_integral / norm_fucktor,
                     "_norm_factor": norm_fucktor}

    # if you're a weenie and don't wanna normalize or edge weight
    e0 = float(min(energies))
    im_coming_out["total_to_Ef"] = {"window": (e0, 0.0),
                          "int_COHP": _integrate_window(energies, cohp, e0, 0.0)}
    return im_coming_out, (energies, cohp)

if __name__ == "__main__":
    edge_vals = {"valence": {}, "conduction": {}, "total_to_Ef": {}}

    for key in hl_gaps:
        path = f"{root}/{key}/"
        if not os.path.exists(f"{path}/tb_input.txt"):
            print(f"{key} not run!")
            run_cogito(directory=path)
            run_cogito_model(dir=path)
        else:
            print(key)

        my_CoTB = CoTB(path)
        my_CoTB.normalize_params()
        my_CoTB.restrict_params(maximum_dist=15, minimum_value=0.00001)

        # choose k-grid getting along shortest lat vec || to 1D chain
        lengths = np.linalg.norm(my_CoTB._a, axis=1)
        base = 24.0
        GRID = tuple(max(2, int(round(base / L))) for L in lengths)
        uni = COGITO_UNIFORM(my_CoTB, GRID)

        vbm_rel, cbm_rel = band_edges_robust(uni, expected_gap=elec_band_gaps[key])

        avail = {}
        for oi in range(uni.num_orbs):
            el = uni.elements[uni.orbatomnum[oi]]
            avail.setdefault(el, set()).add(str(uni.exactorbtype[oi]))

        def _trim(frag):
            out = {}
            for el, want in frag.items():
                have = sorted(set(want) & avail.get(el, set()))
                if have:
                    out[el] = have
            return out

        inorg, org = _trim(INORG), _trim(ORG)

        res, _ = cross_cohp_window(
            uni, vbm_rel, cbm_rel, window=1.0, inorg=inorg, org=org,
            max_dist=3.2, sigma=0.05, spin=0, norm_mode="window",
        )
        for edge in ("valence", "conduction", "total_to_Ef"):
            edge_vals[edge][key] = res[edge]["int_COHP"]
        print(f"{key}: norm={res['_norm_factor']:.3f}  "
              f"val={res['valence']['int_COHP']:+.4f}  "
              f"cond={res['conduction']['int_COHP']:+.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, edge in zip(axes, ("valence", "conduction")):
        xs, ys = [], []
        for key in hl_gaps:
            if key not in edge_vals[edge]:
                continue
            x, y = hl_gaps[key], edge_vals[edge][key]
            xs.append(x); ys.append(y)
            ax.scatter(x, y, color=COLORS[key], marker=MARKERS.get(key, "o"),
                       s=80, edgecolors="black", linewidths=0.8, label=key)
        ax.set_title(f"{edge.capitalize()} band edge")
        ax.set_xlabel("HOMO-LUMO gap (eV)")
        ax.set_ylabel("Integrated cross-COHP (per C-S bond)")
        if len(xs) >= 2:
            r = np.corrcoef(xs, ys)[0, 1]
            ax.text(0.05, 0.95, f"$R$ = {r:.3f}\n$R^2$ = {r**2:.3f}",
                    transform=ax.transAxes, ha="left", va="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(os.path.join(root, "cohp_edges_vs_hlgap.png"), dpi=600, bbox_inches="tight")
    plt.show()