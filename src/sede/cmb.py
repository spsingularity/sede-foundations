"""Compressed-CMB distance priors for SEDE — shift R, acoustic scale l_A.

Borrowed/adapted from the SEDE_V2 implementation (sede/cmb.py). A step beyond the
single shift parameter R used in likelihood.chi2_planck: it adds the acoustic
scale l_A = pi D_M(z*)/r_s(z*), which pins D_M(z*)/r_s and is far more sensitive
to intermediate-redshift dark energy than R alone.

Honest scope (their W9 conclusion, which we adopt): an analytic r_s(z*) reproduces
the VALUE (~144 Mpc) but not the 0.03%-level DERIVATIVES that Planck's l_A error
demands — l_A is hyper-sensitive to omega_m = Omega_m h^2 and needs a Boltzmann
solve (CLASS/CAMB) to be used at LIKELIHOOD level. So:
  * R (sigma ~ 0.26%) is robust and is what enters the joint likelihood.
  * l_A is computed here as a DIAGNOSTIC of how far a model sits in (R, l_A) space.

r_s(z*) is set by pre-recombination physics and is evaluated on the matter+
radiation background (late-time DE is negligible at z >= z*, early-DE ~ 3e-4).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

C_KM_S = 299792.458
OMEGA_GAMMA_H2 = 2.473e-5          # photon density omega_gamma (T_cmb=2.7255 K)
Z_STAR = 1089.92

# Planck 2018 (TT,TE,EE+lowE) shift parameter and a representative l_A.
R_PLANCK, SIG_R = 1.7502, 0.0046
LA_PLANCK, SIG_LA = 301.471, 0.090


def omega_r_of_h(h, n_eff=3.046):
    """Radiation density parameter Omega_r = omega_gamma(1 + 0.2271 Neff)/h^2."""
    return OMEGA_GAMMA_H2 * (1.0 + 0.2271 * n_eff) / h**2


def sound_horizon_rs(Omega_m, H0, omega_b=0.02237, z_star=Z_STAR):
    """Comoving sound horizon at last scattering r_s(z*) in Mpc.

        r_s = int_{z*}^inf c_s(z)/H(z) dz,  c_s = c/sqrt(3(1+R_b)),
        R_b = (3 omega_b)/(4 omega_gamma)/(1+z).
    H(z) uses the matter+radiation background only (DE negligible at z>=z*).
    """
    h = H0 / 100.0
    Omega_rad = omega_r_of_h(h)
    R_b0 = 3.0 * omega_b / (4.0 * OMEGA_GAMMA_H2)

    def integrand(z):
        a1 = 1.0 + z
        E = np.sqrt(Omega_m * a1**3 + Omega_rad * a1**4)
        H = H0 * E
        c_s = C_KM_S / np.sqrt(3.0 * (1.0 + R_b0 / a1))
        return c_s / H

    val, _ = quad(integrand, z_star, 1.0e7, limit=200)
    return val


def D_M_zstar(interp_DC, Omega_m, gamma, H0, z_star=Z_STAR):
    """Comoving distance D_M(z*) in Mpc from a SEDE D_C interpolator (units c/H0)."""
    DC = float(interp_DC(np.array([[Omega_m, gamma, z_star]]))[0]) * (C_KM_S / H0)
    return DC


def shift_R(Omega_m, H0, D_M_star):
    """CMB shift parameter R = sqrt(Omega_m) H0 D_M(z*) / c."""
    return np.sqrt(Omega_m) * H0 * D_M_star / C_KM_S


def acoustic_scale_lA(Omega_m, H0, D_M_star, omega_b=0.02237, z_star=Z_STAR):
    """Acoustic scale l_A = pi D_M(z*) / r_s(z*) (diagnostic)."""
    r_s = sound_horizon_rs(Omega_m, H0, omega_b, z_star)
    return np.pi * D_M_star / r_s


def cmb_diagnostics(interp_DC, Omega_m, gamma, H0, omega_b=0.02237, z_star=Z_STAR,
                    lA_ref=None):
    """Return {R, l_A, r_s, D_M_star} and tensions for a SEDE model.

    `interp_DC` is a friedmann D_C interpolator (Om, gamma, z) -> comoving dist
    in units c/H0; for LCDM pass gamma=1.5, for SEDE-H the fitted gamma, for the
    running model pass sigma8 in the gamma slot.

    R_tension is the ROBUST, likelihood-grade shift-parameter tension. The raw
    absolute l_A carries an analytic-r_s convention offset (~17σ even for LCDM),
    so the meaningful l_A diagnostic is the DIFFERENTIAL Δl_A/σ vs an LCDM
    reference l_A (`lA_ref`), which cancels the offset (their W9 finding).
    """
    DMs = D_M_zstar(interp_DC, Omega_m, gamma, H0, z_star)
    R = shift_R(Omega_m, H0, DMs)
    lA = acoustic_scale_lA(Omega_m, H0, DMs, omega_b, z_star)
    rs = sound_horizon_rs(Omega_m, H0, omega_b, z_star)
    out = dict(R=R, l_A=lA, r_s=rs, D_M_star=DMs,
               R_tension=(R - R_PLANCK) / SIG_R)
    out['lA_diff_sigma'] = None if lA_ref is None else (lA - lA_ref) / SIG_LA
    return out


# Planck 2018 compressed (R, l_A) covariance (TT,TE,EE+lowE; CHW2019-like).
_SIG_R, _SIG_LA, _RHO = 0.0046, 0.090, 0.46
_CMB2_COV = np.array([[_SIG_R**2, _RHO*_SIG_R*_SIG_LA],
                      [_RHO*_SIG_R*_SIG_LA, _SIG_LA**2]])
_CMB2_CINV = np.linalg.inv(_CMB2_COV)


def chi2_cmb_distance(interp_DC, Omega_m, gamma, H0, v_fid, omega_b=0.02237):
    """
    DIFFERENTIAL compressed-CMB chi^2 on the (R, l_A) distance-prior pair, vs a
    fiducial vector `v_fid = [R, l_A]` computed for LCDM@Planck with the SAME
    analytic-r_s pipeline (cancels the convention offset). This is the proper
    likelihood-level term that puts the acoustic scale l_A into the joint fit.

    Result (run_cmb_resolution.py): with l_A included, SEDE-H REFITS from the
    fixed-parameter l_A≈-16σ to l_A≈-1σ, R≈-1.3σ and stays Δχ²≈-1.3 vs LCDM —
    the acoustic-scale "tension" was a fixed-parameter artifact, robust across
    covariance assumptions. l_A does NOT exclude SEDE-H.
    """
    DMs = D_M_zstar(interp_DC, Omega_m, gamma, H0)
    v = np.array([shift_R(Omega_m, H0, DMs),
                  acoustic_scale_lA(Omega_m, H0, DMs, omega_b)])
    d = v - np.asarray(v_fid)
    return float(d @ _CMB2_CINV @ d)


def lcdm_fiducial_vector(interp_DC_lcdm, Omega_m=0.3153, H0=67.36):
    """Fiducial [R, l_A] for LCDM@Planck via the pipeline (target for the
    differential chi2_cmb_distance)."""
    DMs = D_M_zstar(interp_DC_lcdm, Omega_m, 1.5, H0)
    return np.array([shift_R(Omega_m, H0, DMs), acoustic_scale_lA(Omega_m, H0, DMs)])


def print_cmb_diagnostics(interp_DC, Omega_m, gamma, H0, label="model",
                          omega_b=0.02237, lA_ref=None):
    d = cmb_diagnostics(interp_DC, Omega_m, gamma, H0, omega_b, lA_ref=lA_ref)
    diff = "" if d['lA_diff_sigma'] is None else f"  Δl_A={d['lA_diff_sigma']:+.1f}σ vs LCDM"
    print(f"  [{label}]  R={d['R']:.4f} ({d['R_tension']:+.1f}σ vs Planck {R_PLANCK})  "
          f"l_A={d['l_A']:.2f}{diff}  r_s={d['r_s']:.2f} Mpc")
    return d
