"""Dimensional (SI) horizon thermodynamics — the quantitative conjugate identity.

Borrowed/adapted from SEDE_V2's thermodynamics module, then combined with our
EOS-closure theorem (Theorem 5D) to a result neither repo had alone:

    f_sat(0) = Omega_DE / (1 - eps/2)        [conjugate identity, Theorem 1]

With the ΛCDM deceleration eps(0)=1.5 Omega_m this gives 0.897 (SEDE_V2's "Test 1").
But the SEDE-H self-consistent Cai-Kim closure FORCES eps(0)=2 Omega_m (Theorem 5D),
and then

    f_sat(0) = (1 - Omega_m) / (1 - Omega_m) = 1   EXACTLY.

So the f_sat(0)=1 normalisation that Theorem 3 imposed BY HAND is a DERIVED
consequence of the conjugate identity plus the dynamical-horizon closure — it is not
a free choice. This removes the "the conjugate identity is merely definitional"
weak point (W4): the identity + dynamical horizon PREDICT f_sat(0)=1.
"""
from __future__ import annotations

import numpy as np

# SI constants
HBAR = 1.054571817e-34      # J s
K_B = 1.380649e-23          # J/K
C_LIGHT = 2.99792458e8      # m/s
G_NEWTON = 6.67430e-11      # m^3 / (kg s^2)
L_PLANCK_SQ = HBAR * G_NEWTON / C_LIGHT**3   # m^2
MPC_M = 3.0856775814913673e22                # m


def H_si(H0_km_s_mpc):
    """H from km/s/Mpc to 1/s."""
    return H0_km_s_mpc * 1.0e3 / MPC_M


def T_AH(H0_km_s_mpc, epsilon):
    """Cai-Kim apparent-horizon temperature (K): T = (hbar H/2pi k_B)(1-eps/2)."""
    H = H_si(H0_km_s_mpc)
    return HBAR * H / (2.0 * np.pi * K_B) * (1.0 - epsilon / 2.0)


def s_AH(H0_km_s_mpc):
    """Bekenstein-Hawking entropy density s_AH = (3/4) k_B H/(l_P^2 c)  [J/K/m^3]."""
    H = H_si(H0_km_s_mpc)
    return 0.75 * K_B * H / (L_PLANCK_SQ * C_LIGHT)


def rho_crit_energy(H0_km_s_mpc):
    """Critical energy density 3 c^2 H^2/(8 pi G)  [J/m^3]."""
    H = H_si(H0_km_s_mpc)
    return 3.0 * C_LIGHT**2 * H**2 / (8.0 * np.pi * G_NEWTON)


def rho_max(H0_km_s_mpc, epsilon):
    """Bousso-saturated DE density rho_max = T_AH s_AH/k_B = (1-eps/2) rho_crit."""
    return T_AH(H0_km_s_mpc, epsilon) * s_AH(H0_km_s_mpc) / K_B


def conjugate_fsat0(Omega_m, eps0=None):
    """
    f_sat(0) = Omega_DE/(1 - eps0/2) from the conjugate identity. If eps0 is None,
    use the SEDE-H self-consistent closure eps0 = 2 Omega_m (Theorem 5D), which
    returns EXACTLY 1.0. Pass eps0 = 1.5*Omega_m to recover the ΛCDM-eps value 0.897.
    """
    if eps0 is None:
        eps0 = 2.0 * Omega_m
    return (1.0 - Omega_m) / (1.0 - eps0 / 2.0)


def fsat0_is_derived(Omega_m=0.315, tol=1e-9):
    """Check that the SEDE-H closure makes f_sat(0)=1 exactly (derived, not assumed)."""
    return abs(conjugate_fsat0(Omega_m) - 1.0) < tol
