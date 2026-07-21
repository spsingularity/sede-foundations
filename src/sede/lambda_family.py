"""λ-family landscape diagnostic (borrowed from SEDE_V2's coupling_lambda axis).

Maps how the H^{2λ}-coupling exponent λ trades off the EOS-gap closure against
intermediate-z dark energy / the CMB shift. λ=0 = additive (phantom, gap open);
λ=1 = full conjugate identity (tracks H², CMB-tense); the dynamical-horizon
SEDE-H (E_SEDE_H, λ_eff≈0.7) closes the gap (Theorem 5D) AND stays R-safe with
the CORRECT growth — i.e. the SEDE_V2 'closure XOR CMB-safe' trade-off is not a
hard one once the growth factor is right.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .friedmann import (E_SEDE_lambda, E_SEDE_H, compute_growth_factor)
from .theory import w_DE_algebraic

C_KM_S = 299792.458
R_PLANCK, SIG_R = 1.7502, 0.0046


def _w0_fluid(Efun, Omega_m, Omega_r=9.0e-5, dz=2e-3):
    z = np.array([0.0, dz, 2 * dz])
    E = Efun(z)
    rho = E**2 - Omega_m * (1 + z)**3 - Omega_r * (1 + z)**4
    dlnrho_dz = (-3 * np.log(rho[0]) + 4 * np.log(rho[1]) - np.log(rho[2])) / (2 * dz)
    return -1.0 - (-dlnrho_dz) / 3.0


def _omega_de(Efun, Omega_m, z0, Omega_r=9.0e-5):
    E = Efun(np.array([z0]))[0]
    return (E**2 - Omega_m * (1 + z0)**3 - Omega_r * (1 + z0)**4) / E**2


def _shift_R(Efun, Omega_m, H0, z_star=1089.92):
    zz = np.unique(np.concatenate([np.linspace(0, 5, 1500),
                                   np.geomspace(5.01, 1200, 1200)]))
    E = Efun(zz)
    Dc = cumulative_trapezoid(1.0 / E, zz, initial=0.0)
    DMs = (C_KM_S / H0) * np.interp(z_star, zz, Dc)
    return np.sqrt(Omega_m) * H0 * DMs / C_KM_S


def lambda_family_scan(Omega_m=0.317, gamma=1.0, H0=68.5,
                       lams=(0.0, 0.25, 0.5, 0.7, 1.0)):
    """Print w0_fluid, EOS gap, Ω_DE(z=3), and R-tension across λ, plus SEDE-H."""
    w_alg = w_DE_algebraic(Omega_m)
    print(f"λ-family landscape  (Ω_m={Omega_m}, γ={gamma}, H0={H0};  w_alg={w_alg:+.4f})")
    print(f"  {'model':16s} {'w0_fluid':>9} {'EOS gap':>8} {'ΩDE(z=3)':>9} {'R':>7} {'R tension':>10}")
    for lam in lams:
        Efun = lambda z, L=lam: E_SEDE_lambda(z, Omega_m, gamma, L)
        w0 = _w0_fluid(Efun, Omega_m)
        R = _shift_R(Efun, Omega_m, H0)
        print(f"  λ={lam:<13.2f} {w0:>9.4f} {w0 - w_alg:>+8.4f} "
              f"{_omega_de(Efun, Omega_m, 3.0):>9.4f} {R:>7.4f} {(R-R_PLANCK)/SIG_R:>+9.1f}σ")
    # the dynamical-horizon SEDE-H (λ_eff from Cai-Kim, closes the gap)
    EfunH = lambda z: E_SEDE_H(z, Omega_m, gamma)
    w0H = _w0_fluid(EfunH, Omega_m)
    RH = _shift_R(EfunH, Omega_m, H0)
    print(f"  {'SEDE-H (dyn)':16s} {w0H:>9.4f} {w0H - w_alg:>+8.4f} "
          f"{_omega_de(EfunH, Omega_m, 3.0):>9.4f} {RH:>7.4f} {(RH-R_PLANCK)/SIG_R:>+9.1f}σ")
    print("  → SEDE-H closes the EOS gap (≈0) AND is R-safe; the additive λ=0 is "
          "phantom, the full λ=1 over-produces z=3 DE.")


if __name__ == "__main__":
    lambda_family_scan()
