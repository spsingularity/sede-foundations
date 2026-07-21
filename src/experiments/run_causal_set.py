"""
Causal-set test of SEDE's dof-COUNTING premise (§8.5, the counting half).

SEDE's Δ is set by how the horizon dof count scales with radius — area R^{d-2} (Δ=0,
Bekenstein–Hawking) vs spatial volume R^{d-1} (Δ=1, SEDE) — NOT by the quantum state (chaos/SYK
settles that, run_syk_scrambling.py / V60). Causal-set QG contains both scalings; we sprinkle a
2+1D causal set and measure them. FINDING: the canonical horizon-entropy object (causal LINKS
across the horizon) is AREA-law; SEDE's Δ=1 needs the BULK element count instead — a non-default
choice even in the QG framework most hospitable to volume-counting.

Output: console summary, results/causal_set.json, results/fig_causal_set.png.
"""
import json
import os
import numpy as np

from sede import causal_set as cs


def main():
    r = cs.run_experiment(rho=6.0, L=10.0, T=10.0, Rs=(3, 4, 5, 6, 7), n_real=6, seed=0)
    Rs, bulk, horiz = r["Rs"], r["bulk"], r["horizon"]
    print("=" * 74)
    print(f"CAUSAL-SET dof COUNTING (2+1D sprinkle, N={r['N']}, rho={r['rho']})")
    print("=" * 74)
    print("  R           : " + " ".join(f"{R:7.1f}" for R in Rs))
    print("  BULK (disk) : " + " ".join(f"{v:7.0f}" for v in bulk)
          + f"   -> exponent {r['exp_bulk']:.2f}  (expect d-1=2, spatial volume)")
    print("  HORIZON(lnk): " + " ".join(f"{v:7.1f}" for v in horiz)
          + f"   -> exponent {r['exp_horizon']:.2f}  (expect d-2=1, area)")
    print(f"\n  EXPONENT GAP (volume - area) = {r['gap']:.2f}   (expect 1 power)")
    print("  Both countings are NATURAL in the SAME causal set:")
    print("   - canonical horizon entropy = LINKS across the horizon ~ AREA  -> Bekenstein-Hawking (Delta=0)")
    print("   - SEDE's Delta=1 needs s_grav = BULK element count ~ VOLUME    -> a non-default choice")
    print("  In physical 3+1D: bulk ~ R^3 (Delta=1) vs links/area ~ R^2 (Delta=0).")
    print("  => the postulate is a dof-COUNTING claim (volume vs area), not a state/chaos claim.")

    # de Sitter confirmation (conformally flat ⟹ same causal structure ⟹ same fork)
    rd = cs.run_experiment_dS(rho=200.0, L=10.0, eta0=1.0, U=8.0, t0=-2.0,
                              Rs=(3, 4, 5, 6, 7), n_real=6, seed=0)
    print("\n" + "-" * 74)
    print(f"de SITTER confirmation (conformal coords, dS measure, N={rd['N']})")
    print("-" * 74)
    print(f"  BULK ∝ R^{rd['exp_bulk']:.2f} (volume)   HORIZON links ∝ R^{rd['exp_horizon']:.2f} (area)"
          f"   gap={rd['gap']:.2f}")
    print("  => dS is conformally flat, so the causal-set area law (links→Gibbons–Hawking) and the")
    print("     bulk volume law BOTH hold exactly as in flat space. The dS cosmological horizon does")
    print("     NOT prefer volume; the canonical (link) count is area. Same fork, correct geometry.")

    r["de_sitter"] = rd
    r["verdict"] = dict(gap_ge_half=bool(r["gap"] > 0.5),
                        bulk_is_volume=bool(1.7 < r["exp_bulk"] < 2.3),
                        horizon_is_area=bool(0.7 < r["exp_horizon"] < 1.4),
                        dS_same_fork=bool(rd["gap"] > 0.5))
    os.makedirs("results", exist_ok=True)
    with open("results/causal_set.json", "w") as fh:
        json.dump(r, fh, indent=2)
    print("\nwrote results/causal_set.json")
    _make_figure(r)
    return r


def _make_figure(r):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(figure skipped: {e})")
        return
    Rs = np.array(r["Rs"], float)
    fig, ax = plt.subplots(1, 1, figsize=(6.2, 4.6))
    ax.loglog(Rs, r["bulk"], "o-", color="C3",
              label=f"BULK (spatial volume) ∝ R^{r['exp_bulk']:.2f}  → Δ=1 (SEDE)")
    ax.loglog(Rs, r["horizon"], "s-", color="C0",
              label=f"HORIZON links (area) ∝ R^{r['exp_horizon']:.2f}  → Δ=0 (Bekenstein–Hawking)")
    # reference slopes
    ax.loglog(Rs, r["bulk"][0] * (Rs / Rs[0])**2, ":", color="C3", lw=1, label="slope 2 (volume)")
    ax.loglog(Rs, r["horizon"][0] * (Rs / Rs[0])**1, ":", color="C0", lw=1, label="slope 1 (area)")
    ax.set_xlabel("horizon radius R (sprinkling units)")
    ax.set_ylabel("degree-of-freedom count")
    ax.set_title("Causal set: two natural horizon counts differ by one power of R\n"
                 "(the Δ=1-vs-Δ=0 fork — SEDE picks the bulk/volume count)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    os.makedirs("output", exist_ok=True)
    fig.savefig("results/fig_causal_set.png", dpi=130)
    plt.close(fig)
    print("wrote results/fig_causal_set.png")


if __name__ == "__main__":
    main()
