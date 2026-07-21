#!/usr/bin/env python3
"""
Guard: the VOLUME formulation (S_grav ∝ V_AH, no Δ) reproduces Barrow Δ=1 IDENTICALLY.

SEDE's dark sector is stated parameter-free as "the gravitational entropy is the coarse-grained
entanglement entropy in the apparent-horizon volume", S_grav ∝ V_AH ⟺ s_grav = const ⟹ ρ_DE ∝ H·f_sat.
This is algebraically the Barrow maximal-deformation case (A^{1+Δ/2}|_{Δ=1}=A^{3/2}=R³=V, λ=0.5).

This script PROVES the relabeling changes nothing: it solves the volume Friedmann fixed point
ALGEBRAICALLY (independent of the repo's ODE solver),
    E_vol(z) = ½[Ω_DE0 f_sat + √((Ω_DE0 f_sat)² + 4·src)],   src = Ω_m(1+z)³ + Ω_r(1+z)⁴,
and compares E(z), w(z) and the dataset χ² (DESI BAO+Lyα, CMB shift R, fσ8) to
fr.E_SEDE_volume (= E_SEDE_barrow(Δ=1)). Differences are at the numerical-method level only.

Run:  python run_volume_equiv.py
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from sede import friedmann as fr, data_loader as dl

Om, Or, GAM = 0.30, 9.0e-5, 1.4964
Ode0 = 1 - Om - Or
zz = np.concatenate([np.linspace(0, 5, 400), np.geomspace(5.01, 60, 200)])
Dg = fr.compute_growth_factor(zz, Om, Or)


def _fsat(z):
    return (1 - np.exp(-GAM * np.interp(z, zz, Dg) ** 2)) / (1 - np.exp(-GAM))


def E_vol_algebraic(z):
    """Independent closed-form solver of the volume fixed point (no ODE)."""
    z = np.atleast_1d(np.asarray(z, float))
    src = Om * (1 + z) ** 3 + Or * (1 + z) ** 4
    b = Ode0 * _fsat(z)
    return 0.5 * (b + np.sqrt(b ** 2 + 4 * src))


def _Dc(Efun, z):
    g = np.linspace(0, float(np.max(z)) + 0.01, 4000)
    return np.interp(z, g, cumulative_trapezoid(1.0 / Efun(g), g, initial=0.0))


def _chi2_bao(Efun):
    z, types, mean, cov = dl.load_desi_dr2(); Ci = np.linalg.inv(cov)
    X = []
    for zz_, t in zip(z, types):
        DM = _Dc(Efun, np.array([zz_]))[0]; DH = 1.0 / Efun(np.array([zz_]))[0]
        X.append({'DM/rd': DM, 'DH/rd': DH, 'DV/rd': (zz_ * DM ** 2 * DH) ** (1 / 3.)}[t])
    X = np.array(X)
    return mean @ Ci @ mean - (X @ Ci @ mean) ** 2 / (X @ Ci @ X)


def _chi2_R(Efun):
    R = np.sqrt(Om) * _Dc(Efun, np.array([dl.PLANCK_Z_STAR]))[0]
    return ((R - dl.PLANCK_R) / dl.PLANCK_R_ERR) ** 2


def _chi2_fs8(Efun):
    zf, fo, fe = dl.load_fss8()
    Dz, fz = fr.compute_growth_model(zf, Om, lambda z: Efun(np.atleast_1d(z)), Or)
    g = fz * Dz; s8 = np.sum(g * fo / fe ** 2) / np.sum(g * g / fe ** 2)
    return np.sum(((fo - s8 * g) / fe) ** 2)


if __name__ == "__main__":
    print("=" * 76)
    print("VOLUME formulation (S∝V_AH, no Δ)  vs  Barrow Δ=1 — equivalence guard")
    print("=" * 76)

    Evol = lambda z: fr.E_SEDE_volume(np.atleast_1d(z), Om, GAM, Or)   # the model API (no Δ)
    Ebar = lambda z: fr.E_SEDE_barrow(np.atleast_1d(z), Om, GAM, Delta=1.0)

    zt = np.concatenate([np.linspace(0, 10, 500), np.geomspace(10.1, 1100, 200)])
    relE = np.max(np.abs(E_vol_algebraic(zt) - Ebar(zt)) / Ebar(zt))   # independent solver vs ODE
    relAPI = np.max(np.abs(Evol(zt) - Ebar(zt)) / Ebar(zt))            # API alias vs Barrow
    print(f"\n  [E(z)] independent algebraic solver vs Barrow(Δ=1): max rel diff = {relE:.2e}")
    print(f"  [E(z)] E_SEDE_volume API vs E_SEDE_barrow(Δ=1):     max rel diff = {relAPI:.2e}")

    print("\n  dataset χ²:            volume        Barrow(Δ=1)     |diff|")
    ok = relE < 1e-4 and relAPI < 1e-12
    for name, fn in [("DESI BAO (+Lyα)", _chi2_bao), ("CMB shift R", _chi2_R), ("fσ8 (Gold-2018)", _chi2_fs8)]:
        cv, cb = fn(E_vol_algebraic), fn(Ebar)
        print(f"  {name:18s} {cv:12.5f} {cb:14.5f}   {abs(cv - cb):.2e}")
        ok = ok and abs(cv - cb) < 1e-2

    print("\n" + "=" * 76)
    print(f"  [{'PASS' if ok else 'FAIL'}] volume reformulation reproduces Barrow Δ=1 (results unchanged)")
    print("=" * 76)
    import sys
    sys.exit(0 if ok else 1)
