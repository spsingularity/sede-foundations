"""
SEDE joint likelihood for MCMC.

Parameters: theta = (Omega_m, H0, gamma, rd, M_B)
  - Omega_m: matter fraction today
  - H0: Hubble constant [km/s/Mpc]
  - gamma: horizon-saturation coupling (usually fixed to gamma_theory)
  - rd: BAO sound horizon [Mpc]  (free parameter, not fixed)
  - M_B: SN Ia absolute magnitude (nuisance, marginalized or fitted)

Uses pre-tabulated E(z), D_C(z) for O(1ms) evaluations.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import solve_triangular

C_KM_S = 299792.458

# Distance modulus zeropoint (H0 in km/s/Mpc absorbed into M_B)
MU_ZEROPOINT = 25.0  # conventional, absorbed into M_B

# ─── Cosmological distances ────────────────────────────────────────────────

def mu_model(z, H0, Omega_m, gamma, interp_E, interp_DC, M_B):
    """
    Distance modulus mu(z) = 5*log10(D_L / 10pc).
    D_L = (1+z) * D_C  with D_C = (c/H0) * int_0^z dz/E(z).
    """
    z     = np.atleast_1d(z)
    pts   = np.column_stack([np.full_like(z, Omega_m),
                              np.full_like(z, gamma),
                              z])
    DC    = interp_DC(pts) * (C_KM_S / H0)   # Mpc
    DL    = (1 + z) * DC                      # luminosity distance Mpc
    mu    = 5.0 * np.log10(DL) + 25.0
    return mu + M_B


def DM_model(z, H0, Omega_m, gamma, rd, interp_DC):
    """Transverse comoving distance DM(z) = D_C(z) [flat universe]."""
    z   = np.atleast_1d(z)
    pts = np.column_stack([np.full_like(z, Omega_m),
                            np.full_like(z, gamma),
                            z])
    DC  = interp_DC(pts) * (C_KM_S / H0)   # Mpc
    return DC


def DH_model(z, H0, Omega_m, gamma, rd, interp_E):
    """Line-of-sight BAO scale DH(z) = c / H(z)."""
    z   = np.atleast_1d(z)
    pts = np.column_stack([np.full_like(z, Omega_m),
                            np.full_like(z, gamma),
                            z])
    E   = interp_E(pts)
    return C_KM_S / (H0 * E)   # Mpc


def DV_model(z, H0, Omega_m, gamma, rd, interp_E, interp_DC):
    """Spherically-averaged BAO: DV(z) = (z * DM(z)^2 * DH(z))^{1/3} [Mpc]."""
    DM = DM_model(np.atleast_1d(z), H0, Omega_m, gamma, rd, interp_DC)
    DH = DH_model(np.atleast_1d(z), H0, Omega_m, gamma, rd, interp_E)
    return (np.atleast_1d(z) * DM**2 * DH) ** (1.0 / 3.0)


# ─── Individual likelihoods ────────────────────────────────────────────────

def chi2_desi(theta, data, interp_E, interp_DC):
    """DESI DR2 BAO chi^2 with full 13x13 covariance."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    z_arr, types, mean, cov     = data['desi']

    pred = np.zeros(len(z_arr))
    for k, (z, t) in enumerate(zip(z_arr, types)):
        if t == 'DM/rd':
            pred[k] = DM_model(np.array([z]), H0, Omega_m, gamma, rd, interp_DC)[0] / rd
        elif t == 'DH/rd':
            pred[k] = DH_model(np.array([z]), H0, Omega_m, gamma, rd, interp_E)[0] / rd
        elif t == 'DV/rd':
            pred[k] = DV_model(np.array([z]), H0, Omega_m, gamma, rd, interp_E, interp_DC)[0] / rd
        else:
            pred[k] = DM_model(np.array([z]), H0, Omega_m, gamma, rd, interp_DC)[0] / rd

    delta = mean - pred
    try:
        cov_inv = np.linalg.solve(cov, np.eye(len(delta)))
        return float(delta @ cov_inv @ delta)
    except Exception:
        return 1e12


def chi2_pantheon(theta, data, interp_E, interp_DC):
    """Pantheon+ chi^2 with full 1580×1580 STAT+SYS covariance.
    Uses pre-computed Cholesky factor if available (set in load_all_data)."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    z_arr, mu_obs, cov = data['pantheon']
    if z_arr is None:
        return 0.0

    mu_pred = mu_model(z_arr, H0, Omega_m, gamma, interp_E, interp_DC, M_B)
    delta   = mu_obs - mu_pred
    L = data.get('pantheon_L')
    if L is not None:
        v = solve_triangular(L, delta, lower=True, check_finite=False)
        return float(v @ v)
    try:
        L = np.linalg.cholesky(cov)
        v = solve_triangular(L, delta, lower=True, check_finite=False)
        return float(v @ v)
    except Exception:
        return float(delta @ np.linalg.solve(cov, delta))


def chi2_des(theta, data, interp_E, interp_DC):
    """DES-Dovekie chi^2 with analytic M marginalization."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    z_arr, mu_obs, icov = data['des']
    if z_arr is None:
        return 0.0

    mu_pred = mu_model(z_arr, H0, Omega_m, gamma, interp_E, interp_DC, M_B=0.0)
    delta   = mu_obs - mu_pred
    Cinv_d  = icov @ delta
    ones    = np.ones(len(delta))
    Cinv_1  = icov @ ones
    A  = delta @ Cinv_d
    B  = ones  @ Cinv_d
    C  = ones  @ Cinv_1
    return float(A - B**2 / C + np.log(C / (2 * np.pi)))


def chi2_union3(theta, data, interp_E, interp_DC):
    """Union3 chi^2 with anchor marginalization."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    z_arr, mu_obs, icov = data['union3']
    if z_arr is None:
        return 0.0

    mu_pred = mu_model(z_arr, H0, Omega_m, gamma, interp_E, interp_DC, M_B=0.0)
    delta   = mu_obs - mu_pred
    delta   = delta - np.median(delta)
    return float(delta @ icov @ delta)


def chi2_planck(theta, data, interp_E, interp_DC):
    """Planck compressed CMB: H0, Omega_m, shift parameter R."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    H0_pl, sH0, Om_pl, sOm, R_pl, sR, z_star = data['planck']

    # Shift parameter: R = sqrt(Omega_m H0^2) * D_C(z_star) / c
    pts  = np.array([[Omega_m, gamma, z_star]])
    DC_s = interp_DC(pts)[0] * (C_KM_S / H0)
    R_pred = np.sqrt(Omega_m) * H0 * DC_s / C_KM_S

    chi2_H0 = ((H0 - H0_pl) / sH0)**2
    chi2_Om = ((Omega_m - Om_pl) / sOm)**2
    chi2_R  = ((R_pred - R_pl) / sR)**2
    return float(chi2_H0 + chi2_Om + chi2_R)


def chi2_shoes(theta, data, *args):
    """SH0ES H0 prior."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    H0_S, sH0_S = data['shoes']
    return float(((H0 - H0_S) / sH0_S)**2)


def chi2_moresco(theta, data, interp_E, interp_DC):
    """Moresco cosmic chronometers with full covariance."""
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    z_cc, H_obs, cov_cc = data['moresco']

    pts    = np.column_stack([np.full(len(z_cc), Omega_m),
                               np.full(len(z_cc), gamma),
                               z_cc])
    E_arr  = interp_E(pts)
    H_pred = H0 * E_arr
    delta  = H_obs - H_pred
    try:
        cov_inv = np.linalg.solve(cov_cc, np.eye(len(delta)))
        return float(delta @ cov_inv @ delta)
    except Exception:
        return float(np.sum((delta / np.sqrt(np.diag(cov_cc)))**2))


def chi2_fss8(theta, data, interp_E, interp_DC, D_func=None):
    """
    fσ8 RSD growth rate:  fσ8(z) = f(z)·σ8·D(z).

    If data['growth'] = (interp_D, interp_f) is set (the model's SELF-CONSISTENT
    growth, tabulated from compute_growth_model with the model's own E(z)), both
    f(z) and D(z) come from there — exact, model-consistent. Otherwise it falls
    back to the legacy f≈Ω_m(z)^0.55 with ΛCDM-background D(z), which mixed a
    SEDE E(z) into f but a ΛCDM D, an inconsistency that biased fσ8.
    """
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta
    z_f, fs8_obs, fs8_err = data['fss8']

    pts = np.column_stack([np.full(len(z_f), Omega_m),
                           np.full(len(z_f), gamma),
                           z_f])

    growth = data.get('growth')
    if growth is not None and growth[0] is not None:
        interp_D, interp_f = growth
        D_arr = interp_D(pts)
        f_arr = interp_f(pts)
    else:
        # legacy fallback
        if D_func is None:
            from .friedmann import compute_growth_factor
            D_arr = compute_growth_factor(z_f, Omega_m)
        else:
            D_arr = D_func(z_f, Omega_m)
        E_arr = interp_E(pts)
        Om_z  = Omega_m * (1 + z_f)**3 / E_arr**2
        f_arr = Om_z ** 0.55

    fs8_pred = f_arr * sigma8 * D_arr
    delta    = fs8_obs - fs8_pred
    return float(np.sum((delta / fs8_err)**2))


# ─── Joint log-likelihood ────────────────────────────────────────────────────

def log_likelihood(theta, data, interp_E, interp_DC,
                   use_shoes=True, use_pantheon=True,
                   use_des=False, use_union3=False,
                   use_desi=True, use_planck=True,
                   use_moresco=True, use_fss8=True):
    """
    Total log-likelihood = -0.5 * sum(chi^2_sector).
    """
    Omega_m, H0, gamma, rd, M_B, sigma8 = theta

    # Hard priors
    if not (0.1 < Omega_m < 0.6):   return -np.inf
    if not (50 < H0 < 90):           return -np.inf
    if not (0.5 < gamma < 10.0):     return -np.inf   # cap at GAMMA_GRID max
    if not (100 < rd < 200):         return -np.inf
    if not (-20.5 < M_B < -18.0):    return -np.inf
    if not (0.5 < sigma8 < 1.1):     return -np.inf

    total_chi2 = 0.0
    if use_desi     and data.get('desi'):
        total_chi2 += chi2_desi(theta, data, interp_E, interp_DC)
    if use_pantheon and data.get('pantheon') and data['pantheon'][0] is not None:
        total_chi2 += chi2_pantheon(theta, data, interp_E, interp_DC)
    if use_des      and data.get('des') and data['des'][0] is not None:
        total_chi2 += chi2_des(theta, data, interp_E, interp_DC)
    if use_union3   and data.get('union3') and data['union3'][0] is not None:
        total_chi2 += chi2_union3(theta, data, interp_E, interp_DC)
    if use_planck   and data.get('planck'):
        total_chi2 += chi2_planck(theta, data, interp_E, interp_DC)
    if use_shoes    and data.get('shoes'):
        total_chi2 += chi2_shoes(theta, data)
    if use_moresco  and data.get('moresco'):
        total_chi2 += chi2_moresco(theta, data, interp_E, interp_DC)
    if use_fss8     and data.get('fss8'):
        total_chi2 += chi2_fss8(theta, data, interp_E, interp_DC)

    return -0.5 * total_chi2


def log_prior(theta):
    """Flat priors within bounds (already enforced in log_likelihood)."""
    return 0.0


def log_posterior(theta, data, interp_E, interp_DC, **kwargs):
    ll = log_likelihood(theta, data, interp_E, interp_DC, **kwargs)
    if not np.isfinite(ll):
        return -np.inf
    return ll + log_prior(theta)
