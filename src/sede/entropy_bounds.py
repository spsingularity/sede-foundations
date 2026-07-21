"""
Entropy-bound consistency audit of SEDE's volume (Δ=1) dof-counting.

Question (Route C / §8.5): is volume-law horizon entropy S∝R³ even ALLOWED by the established
entropy bounds, or do they forbid it (⟹ answer "area", SEDE wrong)?

We compare, in Planck units (l_P=1, so R is in Planck lengths):
  - holographic / Bekenstein–Hawking bound : S ≤ A/4 = π R²            (smooth 2-area)
  - Bousso covariant entropy bound          : S(light-sheet) ≤ A/4      (same scaling)
  - Bekenstein bound                        : S ≤ 2π R E
  - SEDE volume / Barrow Δ=1                 : S = (A/A0)^{3/2} ∝ R³
against the smooth-area bounds AND against the fractal-horizon reinterpretation (true area of a
space-filling d_H=3 surface ∝ R³, so S = A_fractal/4 saturates the bound w.r.t. the *true* area).

FINDING (honest): vs the SMOOTH-area bounds, S∝R³ overshoots by ∝R (~10⁶¹ at the cosmic horizon)
— so volume-counting is NOT allowed for a smooth horizon. The ONLY consistent escape is that the
horizon is a genuine d_H=3 fractal (true area ∝ R³), which is the postulate restated, not an
independent rescue. ⟹ the bounds neither forbid nor derive volume-counting; they convert it into
"the horizon is space-filling (d_H=3)". Empirics (Δ measurement) remain decisive.
"""
import numpy as np

# cosmic apparent horizon radius in Planck units: R_H = c/(H0) / l_P ≈ 1/ (H0 in Planck units)
R_H_OVER_LP = 8.0e60          # R_H/l_P ≈ M_P/H0 (the seam factor)


def smooth_area(R):
    """Standard geometric area of a 2-sphere of radius R (Planck units)."""
    return 4.0 * np.pi * R**2


def S_holographic(R):
    """Bekenstein–Hawking / holographic bound S = A/4 (smooth area)."""
    return smooth_area(R) / 4.0


def S_volume(R, A0=4.0):
    """SEDE Barrow Δ=1 horizon entropy S = (A/A0)^{3/2} ∝ R³ (A0=4 ⟹ S_BH=(A/4)^1)."""
    return (smooth_area(R) / A0)**1.5


def S_bekenstein(R):
    """Bekenstein bound 2πRE with horizon energy E ~ ρ_DE·V ~ H²·(4/3)πR³, H=1/R (Planck units)."""
    H = 1.0 / R
    rho_DE = H**2                                  # ρ_crit ~ M_P²H² = H² (M_P=1)
    E = rho_DE * (4.0 / 3.0) * np.pi * R**3
    return 2.0 * np.pi * R * E


def fractal_true_area(R, A0=4.0):
    """The 'true' (Hausdorff) area for which S_volume = A_fractal/4 — i.e. the area a space-filling
    (d_H=3) horizon must have. Returns A_fractal and its scaling exponent in R."""
    return 4.0 * S_volume(R, A0)


def audit(R=R_H_OVER_LP):
    """Violation factors of volume-law entropy vs each smooth-area bound, and the fractal escape."""
    Sv = S_volume(R)
    out = {
        "R_over_lP": R,
        "S_volume": Sv,
        "S_holographic": S_holographic(R),
        "S_bekenstein": S_bekenstein(R),
        "viol_holographic": Sv / S_holographic(R),     # ∝ R  (~10⁶¹)
        "viol_bekenstein": Sv / S_bekenstein(R),        # ∝ R
    }
    # scalings (fit exponents over a range of R)
    Rs = np.logspace(1, 6, 30)
    out["exp_S_volume"] = float(np.polyfit(np.log(Rs), np.log(S_volume(Rs)), 1)[0])
    out["exp_S_holo"] = float(np.polyfit(np.log(Rs), np.log(S_holographic(Rs)), 1)[0])
    out["exp_viol"] = float(np.polyfit(np.log(Rs), np.log(S_volume(Rs) / S_holographic(Rs)), 1)[0])
    out["exp_fractal_area"] = float(np.polyfit(np.log(Rs), np.log(fractal_true_area(Rs)), 1)[0])
    # fractal escape: S_volume = A_fractal/4 EXACTLY by construction ⟹ bound respected vs true area
    out["fractal_escape_exact"] = bool(np.allclose(S_volume(Rs), fractal_true_area(Rs) / 4.0))
    return out
