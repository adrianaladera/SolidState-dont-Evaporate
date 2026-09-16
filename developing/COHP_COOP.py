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

root = "/data/NFS/potato/aladera/COGITO/"
elec_band_gaps = {"2,6-dimethyl":1.546, 
          "1naphthyl":1.735, 
          "3methoxy":1.938, 
          "2MMB":1.780, 
          "2methoxy":1.832,
          "2butane":2.175,
          "2propane":2.140,
        #   "gal-hydrated":2.832,
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
        #   "gal-hydrated":5.189,
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
        #   "gal-hydrated":"#00FF38",
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
        #   "gal-hydrated":"^",
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

    n_occ_k  = np.sum(eigs2d <= ef + tol, axis=1) # per-k occupied count at or below ef to some tol
    n_occ    = int(np.median(n_occ_k)) # robust estimate of true occupied cunt
    constant = bool(np.all(n_occ_k == n_occ)) # does each k have the same cunt or naw; F it metal or ef sux

    # split by max occupied set for VBM and min unoccupied cunt set for CBM
    occ, unocc = eigs2d[eigs2d <= ef + tol], eigs2d[eigs2d > ef + tol]
    vbm_abs, cbm_abs = float(occ.max()), float(unocc.min())     # Fermi-split
    # or band counting to get BG w/o ref to ef (just needs occupied count)
    s = np.sort(eigs2d, axis=1)
    vbm_bc, cbm_bc = float(s[:, n_occ - 1].max()), float(s[:, n_occ].min())  # count

    vbm_rel, cbm_rel = vbm_abs - ef, cbm_abs - ef # shift
    gap = cbm_rel - vbm_rel

    print(f"  eigvals range : [{eigs2d.min():.3f}, {eigs2d.max():.3f}] (shifted)")
    print(f"  E_F (eig frame) = {ef:.4f}  [uni.efermi={uni.efermi:.4f} + shift={shift:.4f}]")
    print(f"  n_occ = {n_occ}  (constant across k: {constant})")
    print(f"  gap: fermi-split = {gap:.4f} eV | band-count = {cbm_bc - vbm_bc:.4f} eV")

    # check to see if this shit matches the actual electronic band gap from DFT
    if not constant:
        print("  WARNING: occupied count varies across k — metallic or E_F still off.")
    assert gap > 0.05, (f"Gap {gap:.4f} eV ~ 0 -> E_F in wrong frame; check energy_shift.")
    if expected_gap is not None and abs(gap - expected_gap) > 0.1:
        print(f"  WARNING: gap {gap:.4f} eV != expected {expected_gap} eV.")
    return vbm_rel, cbm_rel

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

def cross_cohp_window(uni, vbm_rel, cbm_rel, window=1.0,
                      inorg=INORG, org=ORG,
                      max_dist=3.2, sigma=0.05, spin=0):
    """Inorganic<->organic pCOHP integrated over the valence and
    conduction edge windows. vbm_rel/cbm_rel are the band edges
    RELATIVE TO E_F (COGITO puts E_F at 0), i.e. vbm_abs - efermi.
    3.2 Å for max dist to filter out distant enough non-interacting
    inorg-org pairs. spin=0 bc these hoes aren't spin-polarized.
    """
    energies, cohp_by_spin = COGITO_UNIFORM.get_COHP_DOS(
        uni, orbs=[inorg, org], max_dist=max_dist, sigma=sigma,
        include_onsite=False, save_plot=False,
    )
    cohp = np.asarray(cohp_by_spin[spin])

    windows = {"valence":    (vbm_rel - window, vbm_rel),
               "conduction": (cbm_rel, cbm_rel + window)}
    out = {}
    for name, (lo, hi) in windows.items():
        out[name] = {"window": (float(lo), float(hi)),
                     "int_COHP": _integrate_window(energies, cohp, lo, hi)}
    # to-E_F total for context (integrate the same curve up to 0)
    e0 = float(min(energies))
    out["total_to_Ef"] = {"window": (e0, 0.0),
                          "int_COHP": _integrate_window(energies, cohp, e0, 0.0)}
    return out, (energies, cohp)

if __name__ == "__main__":
    cohp_vals = {}
    for key in hl_gaps:
        path = f"{root}/{key}/"
        if not os.path.exists(f"{path}/tb_input.txt"):
            print(f"{key} not run!")
            run_cogito(directory=path)
            run_cogito_model(dir=path)
        else:
            print(key)

        # 1. build the TB model from a directory that has run COGITO-
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

        # 3. band edges on the E_F=0 axis (frame-aware + cross-checked)-
        vbm_rel, cbm_rel = band_edges_robust(uni, expected_gap=elec_band_gaps[key])
        print(f"VBM_rel={vbm_rel:.4f}  CBM_rel={cbm_rel:.4f}  gap={cbm_rel-vbm_rel:.4f} eV")

        # 4. trim fragment dicts to orbitals the basis carries
        # search for all AgS (spd) and organics, but not all structures have the same orb types
        # i.e. some projections lack N or Ag(d) or whatever
        avail = {}
        for oi in range(uni.num_orbs):
            el = uni.elements[uni.orbatomnum[oi]]
            avail.setdefault(el, set()).add(str(uni.exactorbtype[oi]))
        print("available orbtypes:", {k: sorted(v) for k, v in avail.items()})

        def _trim(frag):
            out = {}
            for el, want in frag.items():
                have = sorted(set(want) & avail.get(el, set()))
                missing = set(want) - avail.get(el, set())
                if have:    out[el] = have
                if missing: print(f"  dropping {el}:{sorted(missing)} (not in basis)")
            return out

        inorg, org = _trim(INORG), _trim(ORG)
        print(f"INORG used: {inorg}\nORG   used: {org}")

        # 5. the metric: cross-sublattice COHP over the edge windows
        res, _ = cross_cohp_window(
            uni, vbm_rel, cbm_rel, window=1.0, inorg=inorg, org=org,
            max_dist=3.2, sigma=0.05, spin=0,
        )
        # prints valence and conduction int COHP
        for edge, r in res.items():
            print(f"{edge:>12}: int_COHP = {r['int_COHP']:+.4f}  window = {r['window']}")

        # 6. consistency / sign check: curve-to-E_F vs TB-matrix ICOHP-
        ic     = CoTB.get_ICOHP(uni, spin=0)
        idx_in = CoTB.atmorb_dict_to_ind(uni, inorg)
        idx_or = CoTB.atmorb_dict_to_ind(uni, org)
        cross_icohp = float(ic[np.ix_(idx_in, idx_or)].sum().real)

        curve_to_Ef = res["total_to_Ef"]["int_COHP"]
        cohp_vals[key] = curve_to_Ef
        print(f"\ntotal_to_Ef (curve)  = {curve_to_Ef:+.4f}") # the actual ICOHP val
        print(f"cross ICOHP (matrix) = {cross_icohp:+.4f}")
        if abs(cross_icohp) > 1e-9:
            print(f"ratio (curve/matrix) = {curve_to_Ef/cross_icohp:+.3f}  "
                "[+1 consistent | -1 sign flip | ±2 double/half-count]")

    # 7. scatter: integrated cross-COHP vs. HOMO-LUMO gap, per structure--
    fig, ax = plt.subplots()
    for key in hl_gaps:
        ax.scatter(hl_gaps[key], cohp_vals[key],
                   color=COLORS[key], marker=MARKERS[key], label=key,
                   edgecolors="black", linewidths=0.8)
    ax.set_xlabel("HOMO-LUMO gap (eV)")
    ax.set_ylabel("Integrated cross-COHP to $E_F$")
    ax.legend(loc="best", fontsize="small")

    x = np.array([hl_gaps[key] for key in hl_gaps])
    y = np.array([cohp_vals[key] for key in hl_gaps])
    r = np.corrcoef(x, y)[0, 1]
    ax.text(0.05, 0.95, f"$R$ = {r:.3f}\n$R^2$ = {r**2:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize="small",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    fig.tight_layout()
    fig.savefig(os.path.join(root, "cohp_vs_hlgap.png"), dpi=600)

