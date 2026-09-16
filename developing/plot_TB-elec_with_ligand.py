"""plot_TB-elec_with_ligand.py
By Adriana J. Ladera

Combines the three correlation plots produced separately by:
    - inorg_org_mixing-metric.py  (inorganic/organic PDOS mixing ratio, conduction edge)
    - COHP_COOP.py                 (integrated cross-sublattice COHP to E_F)
    - extract_lig_bonding.py       (C-S iCOHP bond-strength descriptor)

into a single 1x3 figure, all plotted against the HOMO-LUMO gap of the
protonated thiol ligand (gas phase). Each panel reuses the COLORS/MARKERS
dicts, iteration scheme, and plotting instructions (scatter + black marker
outline) from its source file.
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


if __name__ == "__main__":
    fig, (ax_mix, ax_cohp, ax_lig) = plt.subplots(1, 3, figsize=(16, 4))

    # =========================================================================
    # panel 1: inorg_org_mixing-metric.py -- conduction-edge mixing ratio
    # =========================================================================
    mix_x, mix_y = [], []
    for key, gap in mixing_mod.hl_gaps.items():
        vasprun_path = os.path.join(ROOT_1D_RAINBOW, key, "dos", "vasprun.xml")
        if not os.path.exists(vasprun_path):
            print(f"skipping {key!r}: no vasprun.xml found at {vasprun_path}")
            continue

        vr = mixing_mod.Vasprun(vasprun_path)
        cdos = vr.complete_dos
        result = mixing_mod.mixing_metric(cdos)

        ratio = result["conduction"]["ratio_in_over_or"]
        mix_x.append(gap)
        mix_y.append(ratio)
        ax_mix.scatter(
            gap,
            ratio,
            color=mixing_mod.COLORS[key],
            marker=mixing_mod.MARKERS[key],
            s=80,
            edgecolors="black",
            linewidths=0.5,
            label=key,
        )

    ax_mix.set_title("Conduction band edge")
    ax_mix.set_xlabel("HOMO-LUMO gap (eV)")
    ax_mix.set_ylabel("I$_{inorganic}$ / I$_{organic}$")
    ax_mix.legend(loc="best", fontsize="small")

    if len(mix_x) >= 2:
        r_mix = np.corrcoef(mix_x, mix_y)[0, 1]
        ax_mix.text(
            0.95, 0.95,
            f"$r$ = {r_mix:.3f}\n$r^2$ = {r_mix**2:.3f}",
            transform=ax_mix.transAxes,
            ha="right", va="top",
        )

    # =========================================================================
    # panel 2: COHP_COOP.py -- integrated cross-sublattice COHP to E_F
    # =========================================================================
    cohp_vals = {}
    for key in cohp_mod.hl_gaps:
        path = f"{ROOT_COGITO}/{key}/"
        if not os.path.exists(f"{path}/tb_input.txt") and not os.path.exists(f"{path}/all_unique_bonds.json"):
            print(f"{key} not run!")
            cohp_mod.run_cogito(directory=path)
            cohp_mod.run_cogito_model(dir=path)
        else:
            print(key)

        # 1. build the TB model from a directory that has run COGITO
        my_CoTB = cohp_mod.CoTB(path)
        my_CoTB.normalize_params()
        my_CoTB.restrict_params(maximum_dist=15, minimum_value=0.00001)

        # 2. k-grid: denser along the 1D chain axis
        lengths = np.linalg.norm(my_CoTB._a, axis=1)
        base = 24.0
        suggested = tuple(max(2, int(round(base / L))) for L in lengths)
        print(f"lattice |a_i| = {np.round(lengths, 3)}  ->  suggested GRID = {suggested}")

        GRID = suggested
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

        # 5. the metric: cross-sublattice COHP over the edge windows
        res, _ = cohp_mod.cross_cohp_window(
            uni, vbm_rel, cbm_rel, window=1.0, inorg=inorg, org=org,
            max_dist=3.2, sigma=0.05, spin=0,
        )
        for edge, r in res.items():
            print(f"{edge:>12}: int_COHP = {r['int_COHP']:+.4f}  window = {r['window']}")

        # 6. consistency / sign check: curve-to-E_F vs TB-matrix ICOHP
        ic = cohp_mod.CoTB.get_ICOHP(uni, spin=0)
        idx_in = cohp_mod.CoTB.atmorb_dict_to_ind(uni, inorg)
        idx_or = cohp_mod.CoTB.atmorb_dict_to_ind(uni, org)
        cross_icohp = float(ic[np.ix_(idx_in, idx_or)].sum().real)

        curve_to_Ef = res["total_to_Ef"]["int_COHP"]
        cohp_vals[key] = curve_to_Ef
        print(f"\ntotal_to_Ef (curve)  = {curve_to_Ef:+.4f}")
        print(f"cross ICOHP (matrix) = {cross_icohp:+.4f}")
        if abs(cross_icohp) > 1e-9:
            print(f"ratio (curve/matrix) = {curve_to_Ef / cross_icohp:+.3f}  "
                  "[+1 consistent | -1 sign flip | ±2 double/half-count]")

    # 7. scatter: integrated cross-COHP vs. HOMO-LUMO gap, per structure
    for key in cohp_mod.hl_gaps:
        ax_cohp.scatter(cohp_mod.hl_gaps[key], cohp_vals[key],
                        color=cohp_mod.COLORS[key], marker=cohp_mod.MARKERS[key], label=key,
                        edgecolors="black", linewidths=0.8)
    ax_cohp.set_xlabel("HOMO-LUMO gap (eV)")
    ax_cohp.set_ylabel("Integrated cross-COHP to $E_F$")
    ax_cohp.legend(loc="best", fontsize="small")

    x_cohp = np.array([cohp_mod.hl_gaps[key] for key in cohp_mod.hl_gaps])
    y_cohp = np.array([cohp_vals[key] for key in cohp_mod.hl_gaps])
    r_cohp = np.corrcoef(x_cohp, y_cohp)[0, 1]
    ax_cohp.text(0.05, 0.95, f"$R$ = {r_cohp:.3f}\n$R^2$ = {r_cohp**2:.3f}",
                transform=ax_cohp.transAxes, ha="left", va="top", fontsize="small",
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
        ax_lig.annotate(lab, (x, y), fontsize=8,
                        textcoords="offset points", xytext=(4, 4))

    if lig_mask.sum() >= 2:
        ax_lig.text(0.05, 0.95, f"$R$ = {r_lig:.3f}\n$R^2$ = {r_lig**2:.3f}",
                    transform=ax_lig.transAxes, ha="left", va="top", fontsize="small",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    sign_note = "neg = stronger bonding" if ligbond_mod.SIGNED else "pos = stronger bonding"
    ax_lig.set_xlabel("HOMO-LUMO gap (eV)")
    ax_lig.set_ylabel(f"C-S iCOHP  [{ligbond_mod.REDUCE}]  (eV)   ({sign_note})")
    ax_lig.set_title("C-S bond strength vs HOMO-LUMO gap")
    ax_lig.axhline(0, linewidth=0.8, color="k")
    ax_lig.grid(True, alpha=0.3)

    # =========================================================================
    # save combined figure
    # =========================================================================
    for ax in (ax_mix, ax_cohp, ax_lig):
        for spine in ax.spines.values():
            spine.set_linewidth(2)

    fig.tight_layout()
    fig.savefig("TB_elec_with_ligand.png", dpi=600, bbox_inches="tight")
    plt.show()
