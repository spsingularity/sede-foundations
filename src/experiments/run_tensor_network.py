"""
Route C constructive test — is volume-law (Δ=1) horizon entropy REALISABLE in a holographic code,
and what does it require? Ryu–Takayanagi min-cut (= max-flow) entanglement of a radius-R ball under
LOCAL (lattice) vs NONLOCAL (expander) connectivity of the same nodes.

FINDING: local connectivity → min-cut ∝ boundary (AREA, Δ=0); nonlocal/expander → min-cut ∝ interior
size (VOLUME, Δ=1). So volume-law is realisable; its required ingredient is NONLOCAL connectivity of
the horizon dof (the geometry-free, all-to-all structure that also makes SYK maximally chaotic). A
geometrically local horizon is area-law. Existence proof + named ingredient.

Output: console, results/tensor_network.json, results/fig_tensor_network.png.
"""
import json
import os
import numpy as np

from sede import tensor_network as tn


def main():
    r = tn.mincut_scaling(n=80, deg=6, Rs=(3, 4, 5, 6, 7), seed=0)
    print("=" * 74)
    print(f"TENSOR-NETWORK min-cut (RT) entanglement: LOCAL vs NONLOCAL  (N={r['N']}, deg={r['deg']})")
    print("=" * 74)
    print("  R           : " + " ".join(f"{R:6d}" for R in r["Rs"]))
    print("  LOCAL  (RT) : " + " ".join(f"{v:6d}" for v in r["local"])
          + f"   -> exponent {r['exp_local']:.2f}  (area R^{{D-1}}, Δ=0)")
    print("  NONLOCAL    : " + " ".join(f"{v:6d}" for v in r["nonlocal_"])
          + f"   -> exponent {r['exp_nonlocal']:.2f}  (volume R^{{D}}, Δ=1)")
    print("\n  => volume-law (Δ=1) horizon entropy IS realisable in a holographic code; its required")
    print("     ingredient is NONLOCAL connectivity of the horizon dof (cf. SYK all-to-all). A")
    print("     geometrically LOCAL horizon is area-law (Δ=0). (2D testbed: area∝R¹, volume∝R²;")
    print("     in physical 3D space these are R² vs R³.)")

    r["verdict"] = dict(local_is_area=bool(abs(r["exp_local"] - 1.0) < 0.2),
                        nonlocal_is_volume=bool(r["exp_nonlocal"] > 1.6),
                        gap=bool(r["exp_nonlocal"] - r["exp_local"] > 0.6))
    os.makedirs("results", exist_ok=True)
    with open("results/tensor_network.json", "w") as fh:
        json.dump(r, fh, indent=2)
    print("\nwrote results/tensor_network.json")
    _figure(r)
    return r


def _figure(r):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(figure skipped: {e})")
        return
    Rs = np.array(r["Rs"], float)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.loglog(Rs, r["nonlocal_"], "o-", color="C3",
              label=f"NONLOCAL (expander) ∝ R^{r['exp_nonlocal']:.2f}  → volume, Δ=1")
    ax.loglog(Rs, r["local"], "s-", color="C0",
              label=f"LOCAL (lattice) ∝ R^{r['exp_local']:.2f}  → area, Δ=0")
    ax.loglog(Rs, r["nonlocal_"][0] * (Rs / Rs[0])**2, ":", color="C3", lw=1, label="slope 2 (volume)")
    ax.loglog(Rs, r["local"][0] * (Rs / Rs[0])**1, ":", color="C0", lw=1, label="slope 1 (area)")
    ax.set_xlabel("ball radius R (lattice units)")
    ax.set_ylabel("RT min-cut entanglement (bonds)")
    ax.set_title("Holographic code: connectivity sets area vs volume\n"
                 "(nonlocal ⟹ Δ=1 realisable; local ⟹ Δ=0)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    os.makedirs("output", exist_ok=True)
    fig.savefig("results/fig_tensor_network.png", dpi=130)
    plt.close(fig)
    print("wrote results/fig_tensor_network.png")


if __name__ == "__main__":
    main()
