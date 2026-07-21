"""Coarse-grained phase-space entropy of a virialised NFW halo, S(M) ~ M^p.

Tests whether the structural entropy weight p (and hence gamma = (1/2)dlnSigma_S/dlnsigma8)
is the gravitational self-ENERGY value 5/3 (gamma=1.5, the original Theorem-4 choice),
the extensive thermal value 1 (gamma~0.27), or the data-preferred ~1.26 (gamma~0.74).

Method (per-particle Gibbs/adiabat entropy of an isotropic NFW halo):
  rho(r)   = NFW, concentration c(M) from Duffy+2008 (z=0),
  sigma_v^2(r) from the isotropic Jeans equation,
  s(r)     = (3/2) ln sigma_v^2(r) - ln rho(r)  + s0   (ln phase-space density per particle),
  S(M)     = (M/m) [ <s>(M) ],  <s> = mass-weighted mean of s(r) over the halo,
  p        = d ln S / d ln M.

The shape part <s>_shape(M) is M-dependent through c(M) and the virial scalings; the
additive Sackur-Tetrode constant s0 (set by the DM particle mass / phase-space cell)
sits in the denominator of p = 1 + (d<s>_shape/dlnM)/(<s>_shape + s0), so p depends on
the DM microphysics. We report the band of p (and gamma) over plausible s0.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad

G = 4.300917270e-6   # kpc (km/s)^2 / Msun
RHO_CRIT0 = 277.5    # Msun/kpc^3 * h^2  (3H0^2/8piG), h^2 absorbed below via h=0.674
H = 0.674
DELTA = 200.0


def concentration(M):
    """Duffy+2008 (z=0, M200crit), c = A (M/Mpivot)^B."""
    return 5.71 * (M / (2e12 / H)) ** (-0.084)


def nfw_halo(M):
    """Return (r_s, rho_s, R_vir, c) for mass M (Msun) at z=0, Delta=200 rho_crit."""
    rho_crit = RHO_CRIT0 * H**2
    R_vir = (3 * M / (4 * np.pi * DELTA * rho_crit)) ** (1.0 / 3.0)   # kpc
    c = concentration(M)
    r_s = R_vir / c
    mu = np.log(1 + c) - c / (1 + c)
    rho_s = M / (4 * np.pi * r_s**3 * mu)
    return r_s, rho_s, R_vir, c


def shape_entropy(M, n=400):
    """Mass-weighted mean of s(r) = (3/2)ln sigma_v^2 - ln rho over the halo (km/s, Msun/kpc^3)."""
    r_s, rho_s, R_vir, c = nfw_halo(M)
    def rho(r): return rho_s / ((r / r_s) * (1 + r / r_s) ** 2)
    def Mr(r):
        x = r / r_s
        return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x / (1 + x))
    # isotropic Jeans: rho sigma^2 (r) = int_r^inf rho(r') G M(r')/r'^2 dr'
    def integrand_jeans(rp): return rho(rp) * G * Mr(rp) / rp**2
    r_grid = np.geomspace(1e-3 * r_s, R_vir, n)
    sig2 = np.empty(n)
    for i, r in enumerate(r_grid):
        val, _ = quad(integrand_jeans, r, 50 * R_vir, limit=100)
        sig2[i] = val / rho(r)
    s_r = 1.5 * np.log(sig2) - np.log(rho(r_grid))     # per-particle log phase-space density
    w = rho(r_grid) * 4 * np.pi * r_grid**2            # dM/dr weight
    return np.trapezoid(s_r * w, r_grid) / np.trapezoid(w, r_grid)


def entropy_slope(M_grid=None):
    """Return (M, <s>_shape, d<s>/dlnM) over the halo mass range."""
    if M_grid is None:
        M_grid = np.geomspace(1e10, 1e16, 13)
    s = np.array([shape_entropy(M) for M in M_grid])
    dsdlnM = np.gradient(s, np.log(M_grid))
    return M_grid, s, dsdlnM


def p_of_s0(s_shape, dsdlnM, s0):
    """Effective weight p = 1 + (d<s>_shape/dlnM)/(<s>_shape + s0) at a given s0."""
    return 1.0 + dsdlnM / (s_shape + s0)


def horizon_deposited_entropy_scaling(M_grid=None):
    """
    Theorem 4C: the entropy a virialising halo adds to the COSMIC HORIZON, ΔS_AH = E_bind/T_AH,
    with T_AH=H/2π the apparent-horizon (reservoir) temperature — mass-independent at a given
    epoch. So ΔS_AH ∝ E_bind ∝ M^{5/3} (slope ≈ 1.64) ⟹ weight p=5/3 ⟹ γ≈1.50. This is the heat
    E_bind deposited into the horizon (where f_sat lives, Thm 1), NOT the halo's own entropy
    E_bind/T_vir (∝M, the Thm-4B p=1). Returns (M, slope_dS_horizon).
    """
    if M_grid is None:
        M_grid = np.geomspace(1e10, 1e16, 13)
    dS = []
    for M in M_grid:
        r_s, rho_s, R_vir, c = nfw_halo(M)
        mu = np.log(1 + c) - c / (1 + c)
        Ebind = 0.5 * G * M**2 / r_s * (1 - 1.0/(1+c)**2 - 2*np.log(1+c)/(1+c)) / mu**2
        dS.append(Ebind)                # / T_AH ∝ H is mass-independent ⟹ slope = E_bind slope
    return M_grid, float(np.polyfit(np.log(M_grid), np.log(np.abs(dS)), 1)[0])


def binding_entropy_scaling(M_grid=None):
    """
    Gravitational BINDING entropy ΔS = |E_bind|/T_vir of an NFW halo vs mass.
    Demonstrates the key point: the binding ENERGY is super-extensive (E_bind ∝ M^1.64
    ≈ M^{5/3}, the origin of the original p=5/3 weight), but the binding ENTROPY is
    EXTENSIVE (ΔS ∝ M^0.98), because T_vir ∝ M^{2/3} divides out the steep energy.
    So no gravitational entropy reaches the super-extensive p>1.2 the data want — the
    M^{5/3} is an ENERGY, not an entropy. Returns (M, slope_Ebind, slope_dS).
    """
    if M_grid is None:
        M_grid = np.geomspace(1e10, 1e16, 13)
    eb, dS = [], []
    for M in M_grid:
        r_s, rho_s, R_vir, c = nfw_halo(M)
        mu = np.log(1 + c) - c / (1 + c)
        Ebind = 0.5 * G * M**2 / r_s * (1 - 1.0/(1+c)**2 - 2*np.log(1+c)/(1+c)) / mu**2
        Tvir = G * M / R_vir            # ~ sigma_v^2
        eb.append(Ebind); dS.append(Ebind / Tvir)
    sl = lambda y: float(np.polyfit(np.log(M_grid), np.log(np.abs(y)), 1)[0])
    return M_grid, sl(eb), sl(dS)
