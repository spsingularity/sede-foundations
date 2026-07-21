"""
Linear perturbations for SEDE-H: gravitational potential growth, the late-time
ISW source, and S8 — the observables a smooth-dark-energy background predicts at
linear order, beyond the background distances used in the MCMC.

Physics (smooth DE — does NOT cluster on sub-horizon scales; it enters only via
the expansion history E(z), so the matter perturbation obeys the standard growth
equation solved in friedmann.compute_growth_model):

  Poisson + matter:  k^2 Phi = -4 pi G a^2 rho_m delta_m  =>  Phi(z) ∝ D(z)/a = D(z)(1+z).
  Hence  d ln Phi / d ln a = d ln D/d ln a - 1 = f(z) - 1.
  In matter domination f≈1 so Phi is constant (NO ISW); when DE takes over f<1,
  Phi DECAYS, sourcing the late-time integrated Sachs-Wolfe effect.

  ISW temperature:  (dT/T)_ISW = 2 ∫ (∂Phi/∂eta) e^{-tau} deta. With deta=-dz/H,
  a=1/(1+z), Phi∝D(1+z), and ∂Phi/∂eta = Phi (f-1) a H, the integrand reduces
  (model-independent constants aside) to a clean scalar SOURCE:

        S_ISW = ∫_0^{z*} D(z) (f(z) - 1) dz        (negative; |S_ISW| = strength)

  The SEDE-H / LCDM ratio of this source is the linear relative ISW amplitude.
  Planck's ISW measurement constrains this at the tens-of-percent level (the
  SEDE V2 bound is <40%).

  S8 = sigma8 * sqrt(Omega_m/0.3): the weak-lensing clustering amplitude.
"""

import numpy as np
from .friedmann import compute_growth_model, E_LCDM, E_SEDE_volume


def potential_growth(z, Omega_m, E_of_z, Omega_r=9.0e-5):
    """
    Return (Phi, f) at redshifts z for background E_of_z. Phi(z) ∝ D(z)(1+z) is
    the (un-normalised) sub-horizon gravitational potential; f=dlnD/dlna.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    D, f = compute_growth_model(z, Omega_m, E_of_z, Omega_r)
    Phi = D * (1.0 + z)
    return Phi, f


def isw_source(Omega_m, E_of_z, z_max=3.0, n=400, Omega_r=9.0e-5):
    """
    Late-time ISW source scalar  S_ISW = ∫_0^{z_max} D(z)(f(z)-1) dz  for the
    given background. Negative (Phi decays); |S_ISW| is the ISW strength.
    """
    z = np.linspace(0.0, z_max, n)
    D, f = compute_growth_model(z, Omega_m, E_of_z, Omega_r)
    return float(np.trapezoid(D * (f - 1.0), z))


def isw_ratio_sedeH(Omega_m, gamma, z_max=3.0, E_run=None):
    """
    Relative linear ISW amplitude  S_ISW(SEDE-H) / S_ISW(LCDM)  at fixed Omega_m.
    >1 means a LARGER ISW signal than LCDM. `E_run`, if given, overrides the
    background with a callable E(z) (e.g. the running model); otherwise the
    CANONICAL background E_SEDE_volume (Barrow Δ=1, λ=0.5) is used — NOT the λ=1
    E_SEDE_H cousin (which would give w0≈−0.86 instead of the canonical ≈−1.0).
    """
    E_sede = E_run if E_run is not None else (lambda zz: E_SEDE_volume(zz, Omega_m, gamma))
    S_sede = isw_source(Omega_m, E_sede, z_max)
    S_lcdm = isw_source(Omega_m, lambda zz: E_LCDM(zz, Omega_m), z_max)
    return S_sede / S_lcdm, S_sede, S_lcdm


def S8(Omega_m, sigma8):
    """Weak-lensing clustering amplitude S8 = sigma8 sqrt(Omega_m/0.3)."""
    return float(sigma8 * np.sqrt(Omega_m / 0.3))


def fsigma8_curve(z, Omega_m, sigma8, E_of_z, Omega_r=9.0e-5):
    """fsigma8(z) = f(z) sigma8 D(z) for the model's self-consistent growth."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    D, f = compute_growth_model(z, Omega_m, E_of_z, Omega_r)
    return f * sigma8 * D


def sedeH_perturbation_report(Omega_m=0.317, gamma=1.0, sigma8=0.72, E_run=None,
                              label="SEDE-H"):
    """
    Print the linear-perturbation observables for a SEDE-H best fit and the LCDM
    reference: potential growth f(0), the relative ISW amplitude, S8, and the
    fsigma8(z) at the RSD pivots. Returns a dict of the numbers.
    """
    E_sede = E_run if E_run is not None else (lambda zz: E_SEDE_volume(zz, Omega_m, gamma))
    # growth today
    D0, f0 = compute_growth_model(np.array([0.0, 0.5, 1.0]), Omega_m, E_sede)
    Dl, fl = compute_growth_model(np.array([0.0, 0.5, 1.0]), Omega_m,
                                  lambda zz: E_LCDM(zz, Omega_m))
    ratio, S_s, S_l = isw_ratio_sedeH(Omega_m, gamma, E_run=E_run)
    s8_curve = fsigma8_curve(np.array([0.15, 0.5, 1.0]), Omega_m, sigma8, E_sede)
    out = dict(f0_sede=float(f0[0]), f0_lcdm=float(fl[0]),
               isw_ratio=float(ratio), S_sede=S_s, S_lcdm=S_l,
               S8=S8(Omega_m, sigma8))
    print(f"  [{label}]  Om={Omega_m:.3f}  sigma8={sigma8:.3f}")
    print(f"    growth rate f(0): SEDE-H={f0[0]:.4f}  LCDM={fl[0]:.4f}  "
          f"(Delta f/f = {(f0[0]/fl[0]-1)*100:+.2f}%)")
    print(f"    ISW source S_ISW: SEDE-H={S_s:.4f}  LCDM={S_l:.4f}  "
          f"ratio={ratio:.3f}  (|dev|={abs(ratio-1)*100:.1f}%, Planck bound <40%)")
    print(f"    S8 = sigma8*sqrt(Om/0.3) = {out['S8']:.4f}  "
          f"(weak-lensing: KiDS~0.76, DES-Y3~0.78, Planck~0.83)")
    print(f"    fsigma8(z=0.15,0.5,1.0) = {np.round(s8_curve,4)}")
    return out
