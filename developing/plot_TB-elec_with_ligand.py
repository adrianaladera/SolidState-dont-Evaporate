"""plot_TB-elec_with_ligand.py
By Adriana J. Ladera

Combines the correlation plots produced separately by:
    - inorg_org_mixing-metric.py  (inorganic/organic PDOS mixing metric,
                                    valence AND conduction band edges)
    - COHP_COOP.py                 (integrated cross-sublattice COHP,
                                    normalized per C-S bond, valence AND
                                    conduction band edges)
    - extract_lig_bonding.py       (C-S iCOHP bond-strength descriptor)

into a single figure, all plotted against the HOMO-LUMO gap of the
protonated thiol ligand (gas phase). Each panel reuses the COLORS/MARKERS
dicts, iteration scheme, and plotting instructions (scatter + black marker
outline) from its source file, and calls into the current (updated) function
signatures of those files rather than duplicating their logic.
"""

import importlib.util
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import COHP_COOP as cohp_mod
import extract_lig_bonding as ligbond_mod

_mixing_spec = importlib.util.spec_from_file_location(
    "inorg_org_mixing_metric",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "inorg_org_mixing-metric.py"),
)
mixing_mod = importlib.util.module_from_spec(_mixing_spec)
_mixing_spec.loader.exec_module(mixing_mod)

# two distinct data roots: 1D_rainbow (DOS/vasprun runs) vs COGITO (TB/COHP runs)
ROOT_1D_RAINBOW = "/data/NFS/potato/aladera/1D_rainbow/"
ROOT_COGITO = "/data/NFS/potato/aladera/COGITO/"

# mirrors the toggles at the top of inorg_org_mixing-metric.py's __main__
MIX_USE_COGITO = True
MIX_METRIC_KEY = "ratio_in_over_or"   # or "ratio_in_over_or" for the legacy metric

# mirrors COHP_COOP.py's cross_cohp_window(norm_mode=...)
COHP_NORM_MODE = "window" # originally cs_bonds but whatev same number

# TAG = "_Em-mod2_elec"
TAG = "_Em-mod2_HL"

def _trim_orbitals(frag, avail):
    """Same trimming logic as the local `_trim` in COHP_COOP.py's __main__."""
    out = {}
    for el, want in frag.items():
        have = sorted(set(want) & avail.get(el, set()))
        missing = set(want) - avail.get(el, set())
        if have:
            out[el] = have
        if missing:
            print(f"  dropping {el}:{sorted(missing)} (not in basis)")
    return out


def _mix_cdos(key):
    """Build the CompleteDos consumed by mixing_metric(), either from a COGITO
    TB model (default) or straight from a VASP DOS run -- same USE_COGITO
    switch as inorg_org_mixing-metric.py's __main__.
    """
    if MIX_USE_COGITO:
        path = os.path.join(ROOT_COGITO, key)
        path += "/"
        # if not os.path.exists(f"{path}/tb_input.txt") and not os.path.exists(f"{path}/all_unique_bonds.json"):
        #     print(f"{key} not run!")
        #     mixing_mod.run_cogito(directory=path)
        #     mixing_mod.run_cogito_model(dir=path)

        my_CoTB = mixing_mod.CoTB(path)
        my_CoTB.normalize_params()
        my_CoTB.restrict_params(maximum_dist=15, minimum_value=0.00001)

        # denser k-grid along the 1D chain axis (shortest lattice vector)
        lengths = np.linalg.norm(my_CoTB._a, axis=1)
        base = 24.0
        suggested = tuple(max(2, int(round(base / L))) for L in lengths)
        print(f"lattice |a_i| = {np.round(lengths, 3)}  ->  suggested GRID = {suggested}")

        uni = mixing_mod.COGITO_UNIFORM(my_CoTB, suggested)
        return uni.get_pymatgen_completedos(uni)

    vasprun_path = os.path.join(ROOT_1D_RAINBOW, key, "dos", "vasprun.xml")
    if not os.path.exists(vasprun_path):
        return None
    return mixing_mod.Vasprun(vasprun_path).complete_dos


if __name__ == "__main__":
    fig = plt.figure(figsize=(11, 5.5))
    gs = fig.add_gridspec(2, 3)
    ax_mix_val = fig.add_subplot(gs[0, 0])
    ax_mix_cond = fig.add_subplot(gs[1, 0])
    ax_cohp_val = fig.add_subplot(gs[0, 1])
    ax_cohp_cond = fig.add_subplot(gs[1, 1])
    ax_lig = fig.add_subplot(gs[:, 2])

    # =========================================================================
    # panel 1: inorg_org_mixing-metric.py -- valence & conduction mixing metric
    # =========================================================================
    mix_edge_axes = {"valence": ax_mix_val, "conduction": ax_mix_cond}
    mix_data = {"valence": {"x": [], "y": []}, "conduction": {"x": [], "y": []}}

    mix_skipped = []
    for key, gap in mixing_mod.hl_gaps.items():
        print(f"hello {key} are you working")
        try:
            cdos = _mix_cdos(key)
            if cdos is None:
                print(f"skipping {key!r}: no DOS source found")
                continue

            result = mixing_mod.mixing_metric(
                cdos, window=1.0, shape="tanh", interface_to="inorganic", normalize="atoms"
            )
        except Exception as e:
            print(f"[skip] {key}: mixing metric failed -- {e!r}")
            mix_skipped.append(key)
            continue

        for edge, ax in mix_edge_axes.items():
            val = result[edge][MIX_METRIC_KEY]
            mix_data[edge]["x"].append(gap)
            mix_data[edge]["y"].append(val)
            ax.scatter(gap, val, color=mixing_mod.COLORS[key], marker=mixing_mod.MARKERS[key],
                       s=80, edgecolors="black", linewidths=0.5, label=key)

    mix_ylabel = ("% organic character at edge" if MIX_METRIC_KEY == "pct_organic"
                  else "I$_{inorganic}$ / I$_{organic}$")
    for edge, ax in mix_edge_axes.items():
        ax.set_title(f"Mixing metric -- {edge} edge", fontsize="small")
        ax.set_xlabel("HOMO-LUMO gap (eV)")
        ax.set_ylabel(mix_ylabel)

        x, y = mix_data[edge]["x"], mix_data[edge]["y"]
        if len(x) >= 2:
            r = np.corrcoef(x, y)[0, 1]
            ax.text(0.95, 0.95, f"$r$ = {r:.3f}\n$r^2$ = {r**2:.3f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize="small")

    # =========================================================================
    # panel 2: COHP_COOP.py -- cross-sublattice COHP, valence & conduction edges
    #          (normalized per C-S bond via cohp_norm_factor)
    # =========================================================================
    cohp_edge_axes = {"valence": ax_cohp_val, "conduction": ax_cohp_cond}
    cohp_edge_vals = {"valence": {}, "conduction": {}}
    cohp_skipped = []

    for key in cohp_mod.hl_gaps:
        path = f"{ROOT_COGITO}/{key}/"

        try:
            # 1. build the TB model from a directory that has run COGITO
            my_CoTB = cohp_mod.CoTB(path)
            my_CoTB.normalize_params()
            my_CoTB.restrict_params(maximum_dist=15, minimum_value=0.00001)

            # 2. k-grid: denser along the 1D chain axis
            lengths = np.linalg.norm(my_CoTB._a, axis=1)
            base = 24.0
            GRID = tuple(max(2, int(round(base / L))) for L in lengths)
            print(f"lattice |a_i| = {np.round(lengths, 3)}  ->  suggested GRID = {GRID}")
            uni = cohp_mod.COGITO_UNIFORM(my_CoTB, GRID)

            # 3. band edges on the E_F=0 axis
            vbm_rel, cbm_rel = cohp_mod.band_edges_robust(
                uni, expected_gap=cohp_mod.elec_band_gaps.get(key)
            )
            print(f"VBM_rel={vbm_rel:.4f}  CBM_rel={cbm_rel:.4f}  gap={cbm_rel - vbm_rel:.4f} eV")

            # 4. trim fragment dicts to orbitals the basis carries
            avail = {}
            for oi in range(uni.num_orbs):
                el = uni.elements[uni.orbatomnum[oi]]
                avail.setdefault(el, set()).add(str(uni.exactorbtype[oi]))
            print("available orbtypes:", {k: sorted(v) for k, v in avail.items()})

            inorg = _trim_orbitals(cohp_mod.INORG, avail)
            org = _trim_orbitals(cohp_mod.ORG, avail)
            print(f"INORG used: {inorg}\nORG   used: {org}")

            # 5. the metric: cross-sublattice COHP over the edge windows, normalized
            #    by the number of C-S contacts so magnitudes are comparable across
            #    structures (COHP_COOP.py's cohp_norm_factor)
            res, _ = cohp_mod.cross_cohp_window(
                uni, vbm_rel, cbm_rel, window=1.0, inorg=inorg, org=org,
                max_dist=3.2, sigma=0.05, spin=0, norm_mode=COHP_NORM_MODE,
                shape="tanh",
            )
            print(f"{key}: norm={res['_norm_factor']:.3f}  "
                  f"val={res['valence']['int_COHP']:+.4f}  "
                  f"cond={res['conduction']['int_COHP']:+.4f}")
        except Exception as e:
            print(f"[skip] {key}: cross-COHP build failed -- {e!r}")
            cohp_skipped.append(key)
            continue

        cohp_edge_vals["valence"][key] = res["valence"]["int_COHP"]
        cohp_edge_vals["conduction"][key] = res["conduction"]["int_COHP"]

    for edge, ax in cohp_edge_axes.items():
        xs, ys = [], []
        for key in cohp_mod.hl_gaps:
            if key not in cohp_edge_vals[edge]:
                continue
            x, y = cohp_mod.hl_gaps[key], cohp_edge_vals[edge][key]
            xs.append(x); ys.append(y)
            ax.scatter(x, y, color=cohp_mod.COLORS[key], marker=cohp_mod.MARKERS.get(key, "o"),
                       s=80, edgecolors="black", linewidths=0.8, label=key)
        ax.set_title(f"Cross-COHP -- {edge} edge", fontsize="small")
        ax.set_xlabel("HOMO-LUMO gap (eV)")
        ax.set_ylabel("Integrated cross-COHP")
        if len(xs) >= 2:
            r = np.corrcoef(xs, ys)[0, 1]
            ax.text(0.05, 0.95, f"$R$ = {r:.3f}\n$R^2$ = {r**2:.3f}",
                    transform=ax.transAxes, ha="left", va="top", fontsize="small",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    # =========================================================================
    # panel 3: extract_lig_bonding.py -- C-S iCOHP bond-strength descriptor
    # =========================================================================
    lig_labels, lig_gaps, lig_cs = [], [], []

    for key, gap in ligbond_mod.hl_gaps.items():
        d = os.path.join(ROOT_COGITO, key)
        path = os.path.join(d, "all_unique_bonds.json")

        if not os.path.isfile(path):
            if ligbond_mod.GENERATE_IF_MISSING:
                path = ligbond_mod.ensure_bond_json(d)
            else:
                print(f"[skip] {key}: {path} not found")
                continue

        print(f"{key}:")
        val = ligbond_mod.cs_bond_strength(
            path, a="C", b="S",
            reduce=ligbond_mod.REDUCE, signed=ligbond_mod.SIGNED, verbose=True,
        )

        lig_labels.append(key)
        lig_gaps.append(gap)
        lig_cs.append(val)
        print(f"   -> C-S ({ligbond_mod.REDUCE}) = {val:+.4f} eV\n")

    lig_gaps = np.asarray(lig_gaps, dtype=float)
    lig_cs = np.asarray(lig_cs, dtype=float)

    lig_mask = ~np.isnan(lig_cs)
    if lig_mask.sum() >= 2:
        r_lig = np.corrcoef(lig_gaps[lig_mask], lig_cs[lig_mask])[0, 1]
        print(f"\nPearson r (gap vs C-S {ligbond_mod.REDUCE}) = {r_lig:+.3f}")

    for x, y, lab in zip(lig_gaps, lig_cs, lig_labels):
        if np.isnan(y):
            continue
        ax_lig.scatter(x, y, color=ligbond_mod.COLORS[lab], marker=ligbond_mod.MARKERS[lab],
                       s=80, edgecolors="black", linewidths=0.8, label=lab)

    if lig_mask.sum() >= 2:
        ax_lig.text(0.05, 0.95, f"$R$ = {r_lig:.3f}\n$R^2$ = {r_lig**2:.3f}",
                    transform=ax_lig.transAxes, ha="left", va="top", fontsize="small",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    sign_note = "neg = stronger bonding" if ligbond_mod.SIGNED else "pos = stronger bonding"
    ax_lig.set_xlabel("HOMO-LUMO gap (eV)")
    ax_lig.set_ylabel(f"C-S iCOHP  [{ligbond_mod.REDUCE}]  (eV)   ({sign_note})")
    ax_lig.set_title("C-S bond strength vs HOMO-LUMO gap", fontsize="small")

    # =========================================================================
    # save combined figure
    # =========================================================================
    for ax in (ax_mix_val, ax_mix_cond, ax_cohp_val, ax_cohp_cond, ax_lig):
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

    
    fig.set_size_inches(12, 4.5, forward=True)
    fig.tight_layout(rect=[0.13, 0, 1, 1])

    handles, labels = ax_mix_val.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.0, 0.5),
               fontsize="x-small")

    fig.savefig(f"TB_elec_with_ligand{TAG}.png", dpi=600, bbox_inches="tight")
    plt.show()

    if mix_skipped:
        print(f"\nmixing panel skipped {len(mix_skipped)} structure(s): {mix_skipped}")
    if cohp_skipped:
        print(f"COHP panel skipped {len(cohp_skipped)} structure(s): {cohp_skipped}")
