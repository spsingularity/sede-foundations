"""
Barrow black-hole thermodynamics — the CROSS-HORIZON sector of SEDE (Tier 3).

SEDE fits a Barrow deformation Δ=1 on the COSMOLOGICAL horizon (Thms 8–9). If Δ is a
universal property of horizons (the quantum-gravity reading, QUANTUM_GRAVITY_IMPLICATIONS.md),
the SAME Δ must appear on BLACK-HOLE horizons. This module gives the BH-sector predictions
that an independent measurement (GW area-law tests, primordial-BH evaporation) can falsify.

Conventions: geometric units G=c=ℏ=k_B=1; masses in solar masses for areas, in Planck
masses for entropy/temperature (M_sun = 1.0e38 M_Pl). Barrow: S=(A/4)^{1+Δ/2} (Δ=0 → A/4,
Bekenstein–Hawking; Δ=1 → maximal fractal, S∝A^{3/2}).  Saridakis 2020; Barrow 2020.
"""
import numpy as np

MSUN_IN_MPL = 1.0e38          # M_sun / M_Planck (order of magnitude)

def kerr_area(M, a_star=0.0):
    """Kerr horizon area in (M_sun)² units of 16π: A = 8πM²[1+√(1−a*²)]."""
    return 8.0 * np.pi * M**2 * (1.0 + np.sqrt(np.maximum(1.0 - a_star**2, 0.0)))

def bh_entropy(M, a_star=0.0, Delta=0.0):
    """Barrow BH entropy S=(A/4)^{1+Δ/2} in Planck units (M given in M_sun)."""
    A = kerr_area(M * MSUN_IN_MPL, a_star)         # area in Planck units
    return (A / 4.0) ** (1.0 + Delta / 2.0)

def bh_temperature_ratio(M, Delta, a_star=0.0):
    """T_Barrow / T_Hawking for a BH. T = (dS/dM)⁻¹; with S∝A^{1+Δ/2}, S∝M^{2+Δ}
    (Schwarzschild) ⇒ T_B/T_H = 1 / [(1+Δ/2)(A/4)^{Δ/2}].  Hugely suppressed for Δ>0."""
    A = kerr_area(M * MSUN_IN_MPL, a_star)
    return 1.0 / ((1.0 + Delta / 2.0) * (A / 4.0) ** (Delta / 2.0))

def entropy_enhancement_log10(M, Delta, a_star=0.0):
    """log10(S_Barrow / S_BH) = (Δ/2) log10(A/4) — the BH-sector signature size."""
    A = kerr_area(M * MSUN_IN_MPL, a_star)
    return (Delta / 2.0) * np.log10(A / 4.0)

def evaporation_time_ratio_log10(M, Delta, a_star=0.0):
    """log10(t_evap,Barrow / t_evap,Hawking). t ∝ ∫ dM/(A T⁴); with T_B suppressed by
    f=(1+Δ/2)(A/4)^{Δ/2}, the rate ∝ T⁴ ⇒ t scales ∝ f⁴ at fixed M (order of magnitude)."""
    A = kerr_area(M * MSUN_IN_MPL, a_star)
    f = (1.0 + Delta / 2.0) * (A / 4.0) ** (Delta / 2.0)
    return 4.0 * np.log10(f)

def area_theorem_margin(m1, m2, Mf, af_star, a1=0.0, a2=0.0):
    """Hawking area theorem: A_f − (A_1+A_2) ≥ 0 (Δ-independent, area is area).
    Returns (A1, A2, Af, ΔA, fractional margin)."""
    A1 = kerr_area(m1, a1); A2 = kerr_area(m2, a2); Af = kerr_area(Mf, af_star)
    dA = Af - (A1 + A2)
    return A1, A2, Af, dA, dA / (A1 + A2)
