"""
Route C audit — is SEDE's volume (Δ=1) dof-counting ALLOWED by the established entropy bounds?

Compares S_volume ∝ R³ against the smooth-area holographic / Bousso / Bekenstein bounds, and against
the fractal-horizon reinterpretation. Output: console, results/entropy_bounds.json, results/fig_entropy_bounds.png.

FINDING: vs SMOOTH-area bounds, S∝R³ overshoots by ∝R (~10⁶¹ at the cosmic horizon) — volume-counting
is NOT allowed for a smooth horizon. The ONLY escape is a genuine d_H=3 fractal horizon (true area
∝ R³), which is the postulate restated. So the bounds neither forbid nor derive volume-counting; they
convert it into "the horizon is space-filling." Empirics (Δ) stay decisive.
"""
import json
import os
import numpy as np

from sede import entropy_bounds as eb


def main():
    a = eb.audit()
    print("=" * 74)
    print("ENTROPY-BOUND AUDIT of SEDE volume-law (Δ=1) at the cosmic horizon")
    print("=" * 74)
    print(f"  R/l_P (seam factor)        : {a['R_over_lP']:.1e}")
    print(f"  S_volume (Barrow Δ=1)      : {a['S_volume']:.2e}   (∝ R^{a['exp_S_volume']:.2f})")
    print(f"  S_holographic = A/4        : {a['S_holographic']:.2e}   (∝ R^{a['exp_S_holo']:.2f})")
    print(f"  violation vs holographic   : {a['viol_holographic']:.2e}   (∝ R^{a['exp_viol']:.2f} = R/l_P)")
    print(f"  violation vs Bekenstein    : {a['viol_bekenstein']:.2e}")
    print(f"\n  fractal escape: S_volume = A_fractal/4 with A_fractal ∝ R^{a['exp_fractal_area']:.2f}"
          f"  (exact: {a['fractal_escape_exact']})")
    print("  => vs SMOOTH-area bounds, volume-law overshoots by ~R/l_P (~10^61): NOT allowed for a")
    print("     smooth horizon. RESCUED only if the horizon is a genuine d_H=3 fractal (true area ∝ R³)")
    print("     — which IS the postulate. Bounds neither forbid nor derive Δ=1; they translate it to")
    print("     'the horizon is space-filling'. Empirical Δ measurement remains decisive.")

    a["verdict"] = dict(
        violates_smooth_bound=bool(a["viol_holographic"] > 1e3),
        violation_is_linear_in_R=bool(abs(a["exp_viol"] - 1.0) < 0.05),
        fractal_escape_consistent=bool(a["fractal_escape_exact"]),
    )
    os.makedirs("results", exist_ok=True)
    with open("results/entropy_bounds.json", "w") as fh:
        json.dump({k: (v if not isinstance(v, np.generic) else v.item()) for k, v in a.items()}, fh, indent=2)
    print("\nwrote results/entropy_bounds.json")
    _figure(a)
    return a


def _figure(a):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(figure skipped: {e})")
        return
    R = np.logspace(1, 6, 40)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.loglog(R, eb.S_volume(R), "-", color="C3", label="S_volume ∝ R³ (Barrow Δ=1, SEDE)")
    ax.loglog(R, eb.S_holographic(R), "-", color="C0", label="holographic/Bousso bound A/4 ∝ R²")
    ax.loglog(R, eb.S_bekenstein(R), "--", color="C2", label="Bekenstein bound 2πRE ∝ R²")
    ax.loglog(R, eb.fractal_true_area(R) / 4, ":", color="k", lw=2,
              label="A_fractal/4 (d_H=3 true area) = S_volume (escape)")
    ax.set_xlabel("horizon radius R / l_P")
    ax.set_ylabel("entropy")
    ax.set_title("Volume-law entropy vs the area bounds\n(overshoots smooth A/4 by ∝R; saturates the fractal area)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    os.makedirs("output", exist_ok=True)
    fig.savefig("results/fig_entropy_bounds.png", dpi=130)
    plt.close(fig)
    print("wrote results/fig_entropy_bounds.png")


if __name__ == "__main__":
    main()
