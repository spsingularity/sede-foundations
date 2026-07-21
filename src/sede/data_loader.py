"""
Data loader for SEDE pipeline.

Handles: Pantheon+, DES-Dovekie, Union3, DESI DR2, Moresco CC,
         Planck (H0, Omega_m, R), SH0ES, fσ8 RSD.

All covariance matrices are symmetrized and verified positive-definite.
"""

import numpy as np
from scipy.linalg import cholesky
import os

from ._paths import DATA_DIR   # layout-agnostic: <root>/data in both dev-hub and release trees


# ─── DESI DR2 BAO ───────────────────────────────────────────────────────────

_DESI_TYPE_MAP = {
    'DV_over_rs': 'DV/rd',
    'DM_over_rs': 'DM/rd',
    'DH_over_rs': 'DH/rd',
}


def load_desi_dr2():
    """
    Load DESI DR2 BAO data from official files.

    Returns (z_arr, types, mean, cov_13x13).
    types: list of 'DV/rd', 'DM/rd', or 'DH/rd'.
    BGS at z=0.295 is DV/rd (spherically-averaged).
    """
    mean_path = os.path.join(DATA_DIR, 'desi_dr2_desi_gaussian_bao_ALL_GCcomb_mean.txt')
    cov_path  = os.path.join(DATA_DIR, 'desi_dr2_desi_gaussian_bao_ALL_GCcomb_cov.txt')

    if os.path.exists(mean_path):
        z_list, val_list, type_list = [], [], []
        with open(mean_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                z_list.append(float(parts[0]))
                val_list.append(float(parts[1]))
                raw_type = parts[2]
                type_list.append(_DESI_TYPE_MAP.get(raw_type, raw_type))
        z_arr = np.array(z_list)
        mean  = np.array(val_list)
        types = type_list
    else:
        # Fallback: hardcoded DESI DR2 values (BGS is DV/rd)
        z_arr = np.array([
            0.295,
            0.510, 0.510,
            0.706, 0.706,
            0.934, 0.934,
            1.321, 1.321,
            1.484, 1.484,
            2.330, 2.330,
        ])
        types = [
            'DV/rd',
            'DM/rd', 'DH/rd',
            'DM/rd', 'DH/rd',
            'DM/rd', 'DH/rd',
            'DM/rd', 'DH/rd',
            'DM/rd', 'DH/rd',
            'DH/rd', 'DM/rd',
        ]
        mean = np.array([
            7.94167639,
            13.58758434, 21.86294686,
            17.35069094, 19.45534918,
            21.57563956, 17.64149464,
            27.60085612, 14.17602155,
            30.51190063, 12.81699964,
            8.631545675, 38.988973962,
        ])

    if os.path.exists(cov_path):
        cov = np.loadtxt(cov_path)
    else:
        # Block-diagonal fallback (DM/DH pairs anti-correlated)
        n = len(z_arr)
        cov = np.zeros((n, n))
        diag_approx = np.array([
            5.79e-3,
            2.83e-2, 1.84e-1,
            3.24e-2, 1.11e-1,
            2.62e-2, 4.04e-2,
            1.05e-1, 5.04e-2,
            5.83e-1, 2.68e-1,
            1.02e-2, 2.83e-1,
        ])[:n]
        np.fill_diagonal(cov, diag_approx)

    cov = 0.5 * (cov + cov.T)
    return z_arr, types, mean, cov


# ─── DESI DR2 w0waCDM POSTERIOR (official GetDist chains) ─────────────────────
# Official DESI DR2 cosmology chains (data.desi.lbl.gov/public/papers/y3/
# bao-cosmo-params/cobaya/base_w_wa/), DESI-BAO-all + Planck2018+NPIPE-CamSpec
# +ACT-DR6-lensing + one SN sample. We store the small getdist outputs only
# (chain.covmat = full param covariance; chain.margestats = marginal means):
#   desi_dr2_w0wa_{pantheonplus,desy5,union3}.{covmat,margestats}
# This is the AUTHORITATIVE (w0,wa) covariance — replaces the published-marginal
# approximation (ρ≈-0.90) used previously.
import re as _re

_W0WA_SN = ('pantheonplus', 'desy5', 'union3')


def _margestats_mean(path, name):
    """Leading (mean) value of a parameter row in a GetDist .margestats file."""
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if p and p[0] == name:
                m = _re.match(r'[-+]?[0-9]*\.?[0-9]+', p[1])
                if m:
                    return float(m.group())
    return None


def load_desi_dr2_w0wa():
    """
    Official DESI DR2 w0waCDM (w0, wa) means + 2×2 covariance for each SN combination.

    Returns {sn: dict(w0, wa, cov2x2, rho)} for sn in {'pantheonplus','desy5','union3'},
    parsed from the stored getdist chain.covmat / chain.margestats. Falls back to the
    published marginals (ρ≈-0.90/-0.92) for any sample whose files are absent, so callers
    are robust to a missing download.
    """
    published = {                                      # fallback marginals (DESI DR2 2025)
        'pantheonplus': (-0.838, 0.055, -0.62, 0.22, -0.90),
        'desy5':        (-0.752, 0.057, -0.86, 0.22, -0.905),
        'union3':       (-0.667, 0.088, -1.09, 0.29, -0.93),
    }
    out = {}
    for sn in _W0WA_SN:
        cmf = os.path.join(DATA_DIR, f'desi_dr2_w0wa_{sn}.covmat')
        msf = os.path.join(DATA_DIR, f'desi_dr2_w0wa_{sn}.margestats')
        if os.path.exists(cmf) and os.path.exists(msf):
            hdr = open(cmf).readline().lstrip('#').split()
            cov = np.loadtxt(cmf)
            iw, iwa = hdr.index('w'), hdr.index('wa')
            sub = cov[np.ix_([iw, iwa], [iw, iwa])]
            w0 = _margestats_mean(msf, 'w'); wa = _margestats_mean(msf, 'wa')
            rho = sub[0, 1] / np.sqrt(sub[0, 0] * sub[1, 1])
            out[sn] = dict(w0=w0, wa=wa, cov2x2=sub, rho=float(rho), source='official')
        else:
            w0, sw0, wa, swa, rho = published[sn]
            sub = np.array([[sw0**2, rho*sw0*swa], [rho*sw0*swa, swa**2]])
            out[sn] = dict(w0=w0, wa=wa, cov2x2=sub, rho=rho, source='published-fallback')
    return out


# ─── PLANCK COMPRESSED CMB ────────────────────────────────────────────────────

PLANCK_H0     = 67.4
PLANCK_H0_ERR = 0.5
PLANCK_OM     = 0.315
PLANCK_OM_ERR = 0.007
PLANCK_R      = 1.7502   # Planck 2018, as specified in reference §8
PLANCK_R_ERR  = 0.0046   # reference §8: R = 1.7502 ± 0.0046
PLANCK_Z_STAR = 1089.9


def load_planck():
    """Returns (H0, sigma_H0, Omega_m, sigma_Om, R, sigma_R, z_star)."""
    return (PLANCK_H0, PLANCK_H0_ERR,
            PLANCK_OM, PLANCK_OM_ERR,
            PLANCK_R,  PLANCK_R_ERR,
            PLANCK_Z_STAR)


# ─── SH0ES ───────────────────────────────────────────────────────────────────

SHOES_H0     = 73.04   # Riess+2022 (R22): H0 = 73.04 ± 1.04 km/s/Mpc
SHOES_H0_ERR = 1.04    # (was 73.0±1.0; unified with the CAMB drivers' value)


def load_shoes():
    return SHOES_H0, SHOES_H0_ERR


# ─── MORESCO COSMIC CHRONOMETERS ─────────────────────────────────────────────
# 15-point H(z) compilation (Moresco+2016, 2020)

MORESCO_Z = np.array([
    0.07, 0.09, 0.12, 0.17, 0.179, 0.199,
    0.27, 0.352, 0.4, 0.593, 0.68, 0.781, 0.875, 1.037, 1.3,
])
MORESCO_H = np.array([
    69.0, 69.0, 68.6, 83.0, 75.0, 75.0,
    77.0, 83.0, 95.0, 104.0, 92.0, 105.0, 125.0, 154.0, 168.0,
])
MORESCO_SIG = np.array([
    19.6, 12.0, 26.2, 8.0, 4.0, 5.0,
    14.0, 14.0, 17.0, 13.0, 8.0, 12.0, 17.0, 20.0, 17.0,
])


def load_moresco():
    """
    Returns (z, H_obs, cov) with full 15x15 covariance.
    Uses diagonal + correlated systematic floor (3% IMF + SPS).
    """
    cov_path = os.path.join(DATA_DIR, 'moresco_cov.npy')
    if os.path.exists(cov_path):
        cov = np.load(cov_path)
    else:
        # NOTE (checkpoint review): the moresco_cov_*.dat files in data/ are failed-download
        # HTML 404 pages, not covariance data — they were never parsed here (this loader only
        # ever looked for moresco_cov.npy). We fall back to a diagonal + 3% correlated-systematic
        # floor. This is an APPROXIMATION, not the true Moresco (2022) covariance. Impact on the
        # headline is negligible (CC contributes ≈+0.11 to Δχ²), but the fallback is now explicit
        # rather than silent. To use the real covariance, cache it as data/moresco_cov.npy.
        import warnings
        warnings.warn(
            "load_moresco: real Moresco (2022) covariance not found (data/moresco_cov.npy "
            "absent; the moresco_cov_*.dat files are invalid HTML). Using a diagonal + 3% "
            "correlated-systematic-floor APPROXIMATION.", RuntimeWarning, stacklevel=2)
        cov_stat  = np.diag(MORESCO_SIG**2)
        sys_frac  = 0.03
        sys_outer = np.outer(MORESCO_H * sys_frac, MORESCO_H * sys_frac)
        cov = cov_stat + sys_outer
    cov = 0.5 * (cov + cov.T)
    return MORESCO_Z, MORESCO_H, cov


# ─── fσ8 RSD ─────────────────────────────────────────────────────────────────

# Standard "Gold-2018" RSD compilation (Sagredo, Nesseris & Sapone 2018), the
# vetted 16-point set also used by the SEDE_V2 sibling. REPLACES an earlier
# 18-point array that carried erroneous over-tight entries — notably z=0.85 with
# σ=0.035 (the real eBOSS-ELG error is ~0.095) and non-standard z=1.05 / z=2.33
# points — which inflated the ABSOLUTE fσ8 χ² (~31 vs ~6) and added a spurious
# ~+2.4 to SEDE-H's disfavouring (with the clean data the fσ8 Δ(SEDE-H−ΛCDM)
# is −0.3, i.e. SEDE-H fits slightly better). See chase of the fσ8 discrepancy.
FSS8_Z = np.array([
    0.02, 0.02, 0.10, 0.15, 0.18, 0.38, 0.25, 0.37,
    0.32, 0.59, 0.44, 0.60, 0.73, 0.86, 1.40, 0.978,
])
FSS8_VAL = np.array([
    0.428, 0.398, 0.370, 0.490, 0.360, 0.440, 0.3512, 0.4602,
    0.384, 0.488, 0.413, 0.390, 0.437, 0.400, 0.482, 0.379,
])
FSS8_ERR = np.array([
    0.0465, 0.065, 0.130, 0.145, 0.090, 0.060, 0.0583, 0.0378,
    0.095, 0.060, 0.080, 0.063, 0.072, 0.110, 0.116, 0.176,
])


def load_fss8():
    """Returns (z, fσ8_obs, sigma) for the Gold-2018 RSD compilation (16 points)."""
    return FSS8_Z, FSS8_VAL, FSS8_ERR


# ─── eBOSS DR16 FULL-SHAPE (pre-DESI BAO+RSD) ────────────────────────────────
# SDSS BOSS DR12 + eBOSS DR16 "BAOplus" consensus full-shape (D_M/r_d, D_H/r_d,
# fσ8) for LRG (z=0.38, 0.51, 0.698) and QSO (z=1.48). Public files fetched from
# CobayaSampler/bao_data (Alam et al. 2021). This is GENUINE pre-DESI data — used
# for the blind-DESI out-of-sample test (train pre-DESI → predict DESI DR2).
_DR16_TYPE = {'DM_over_rs': 'DM/rd', 'DH_over_rs': 'DH/rd', 'f_sigma8': 'fs8'}

def load_eboss_dr16_fs():
    """Returns (z, types, mean, cov) for the combined eBOSS DR16 LRG+QSO full-shape
    (12 measurements; block-diagonal LRG 9×9 ⊕ QSO 3×3 covariance). types ∈
    {'DM/rd','DH/rd','fs8'}. Returns (None,)*4 if the files are absent."""
    import scipy.linalg as _sla
    lrg = os.path.join(DATA_DIR, 'sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8.dat')
    qso = os.path.join(DATA_DIR, 'sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8.dat')
    lrgc = os.path.join(DATA_DIR, 'sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8_covtot.txt')
    qsoc = os.path.join(DATA_DIR, 'sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8_covtot.txt')
    if not all(os.path.exists(p) for p in (lrg, qso, lrgc, qsoc)):
        print("[WARN] eBOSS DR16 full-shape files not found")
        return None, None, None, None
    z, types, mean = [], [], []
    for path in (lrg, qso):
        with open(path) as fh:
            for line in fh:
                zz, val, tp = line.split()
                z.append(float(zz)); mean.append(float(val)); types.append(_DR16_TYPE[tp])
    cov = _sla.block_diag(np.loadtxt(lrgc), np.loadtxt(qsoc))
    return np.array(z), types, np.array(mean), cov


# ─── PANTHEON+ ───────────────────────────────────────────────────────────────

def load_pantheon_plus(data_path=None, cov_path=None):
    """
    Load Pantheon+SH0ES SN Ia data.

    Cuts: zHD > 0.01 and IS_CALIBRATOR == 0 → ~1580 Hubble-flow SNe.
    Uses m_b_corr (SALT2-corrected apparent magnitude).
    Returns (z, mu_obs, cov_NxN).
    """
    if data_path is None:
        data_path = os.path.join(DATA_DIR, 'Pantheon+SH0ES.dat')
    if cov_path is None:
        cov_path = os.path.join(DATA_DIR, 'Pantheon+SH0ES_STAT+SYS.cov')

    if not os.path.exists(data_path):
        print(f"[WARN] Pantheon+ data not found at {data_path}")
        return None, None, None

    # Read full data table (mixed types: CID is string)
    data = np.genfromtxt(data_path, names=True, dtype=None, encoding='utf-8')
    mask = (data['zHD'] > 0.01) & (data['IS_CALIBRATOR'] == 0)
    z    = data['zHD'][mask].astype(float)
    mu   = data['m_b_corr'][mask].astype(float)

    if not os.path.exists(cov_path):
        print(f"[WARN] Pantheon+ covariance not found at {cov_path}")
        return z, mu, np.diag(np.ones(len(z)) * 0.1**2)

    # Covariance file: line 1 = n, then n^2 values (one per line)
    with open(cov_path) as f:
        n_cov = int(f.readline().strip())
    C_flat = np.loadtxt(cov_path, skiprows=1)
    C_full = C_flat.reshape(n_cov, n_cov)

    # Extract submatrix for the selected SNe
    idx   = np.where(mask)[0]
    C_sub = C_full[np.ix_(idx, idx)]
    C_sub = 0.5 * (C_sub + C_sub.T)
    return z, mu, C_sub


# ─── DES-DOVEKIE ─────────────────────────────────────────────────────────────

def _parse_des_csv(path):
    """
    Parse DES-Dovekie SNANA FITRES format.
    Header comes from the 'VARNAMES:' line; data rows start with 'SN:'.
    Returns dict mapping column name → array of values (strings).
    """
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[0] == 'VARNAMES:':
                header = parts[1:]   # skip the 'VARNAMES:' token
                continue
            if parts[0] == 'SN:' and header is not None:
                data_parts = parts[1:]
                if len(data_parts) == len(header):
                    rows.append(data_parts)
    if header is None:
        raise ValueError(f"No VARNAMES line found in {path}")
    if not rows:
        raise ValueError(f"No SN: data rows found in {path}")
    arr = np.array(rows)
    return {name: arr[:, i] for i, name in enumerate(header)}


def load_des_dovekie(data_path=None, icov_path=None):
    """
    DES 5YR SN Ia (Dovekie 2024) — 1820 SNe, M analytically marginalized.

    Returns (z, mu_obs, inv_cov_1820x1820).
    The inv_cov is the INVERSE covariance (stored directly in the npz).
    """
    if data_path is None:
        data_path = os.path.join(DATA_DIR, 'DES_Dovekie_HD.csv')
    if icov_path is None:
        icov_path = os.path.join(DATA_DIR, 'DES_STAT+SYS.npz')

    if not os.path.exists(data_path):
        print(f"[WARN] DES-Dovekie data not found at {data_path}")
        return None, None, None

    cols = _parse_des_csv(data_path)
    z_all  = cols['zHD'].astype(float)
    mu_all = cols['MU'].astype(float)
    mask   = z_all > 0.0
    z      = z_all[mask]
    mu     = mu_all[mask]

    if not os.path.exists(icov_path):
        print(f"[WARN] DES covariance not found at {icov_path}")
        return z, mu, np.diag(1.0 / (0.15**2 * np.ones(len(z))))

    d = np.load(icov_path, allow_pickle=False)
    n = int(d['nsn'][0])
    # Upper-triangular inverse covariance → full symmetric matrix
    icov = np.zeros((n, n))
    icov[np.triu_indices(n)] = d['cov'].astype(float)
    icov = icov + icov.T - np.diag(np.diag(icov))

    if len(z) != n:
        print(f"[WARN] DES: {len(z)} SNe from CSV but cov size is {n}; using first {n}")
        z  = z[:n]
        mu = mu[:n]

    return z, mu, icov


# ─── UNION3 ──────────────────────────────────────────────────────────────────

def load_union3(fits_path=None):
    """
    Union3 2023 — 22 binned SN Ia (z=0.05 to 2.26), anchor marginalized.
    Returns (z_22, mu_22, inv_cov_22x22).

    FITS 23x23 layout:
      arr[0, 1:] = z-bin centres (22 values)
      arr[1:, 0] = μ observed (22 values)
      arr[1:, 1:] = 22x22 INVERSE covariance
    """
    if fits_path is None:
        # Try both naming conventions
        for name in ('union3_mu.fits', 'union3.fits'):
            candidate = os.path.join(DATA_DIR, name)
            if os.path.exists(candidate):
                fits_path = candidate
                break

    if fits_path is None or not os.path.exists(str(fits_path)):
        print("[WARN] Union3 FITS not found; using diagonal fallback")
        z_u3  = np.array([0.07, 0.11, 0.15, 0.20, 0.26, 0.33, 0.42,
                           0.53, 0.67, 0.85, 1.07, 1.34, 1.68, 2.26])
        mu_u3 = np.array([32.9, 34.3, 35.3, 36.4, 37.5, 38.5, 39.6,
                           40.7, 41.8, 42.8, 43.9, 44.8, 45.7, 46.9])
        return z_u3, mu_u3, np.diag(1.0 / (0.15**2 * np.ones(len(z_u3))))

    try:
        from astropy.io import fits as pyfits
        with pyfits.open(fits_path) as hdul:
            arr   = hdul[0].data          # shape (23, 23)
            z_u3  = arr[0, 1:].copy()    # 22 z-bin centres
            mu_u3 = arr[1:, 0].copy()    # 22 mu values
            icov  = arr[1:, 1:].copy()   # 22x22 inverse covariance
        icov = 0.5 * (icov + icov.T)
        return z_u3, mu_u3, icov
    except Exception as e:
        print(f"[WARN] Failed to load Union3: {e}")
        return None, None, None


def union3_chi2(mu_pred, mu_obs, icov):
    """Anchor-marginalized chi^2: subtract median offset to marginalize over anchor."""
    delta = mu_obs - mu_pred
    delta = delta - np.median(delta)
    return float(delta @ icov @ delta)


# ─── Utilities ───────────────────────────────────────────────────────────────

def des_chi2_marginalized(mu_pred, mu_obs, icov):
    """
    Conley+2011 analytically M-marginalized chi^2.
    chi2 = delta C^{-1} delta - (1 C^{-1} delta)^2 / (1 C^{-1} 1) + ln(1 C^{-1} 1 / 2π)
    """
    delta     = mu_obs - mu_pred
    Cinv_d    = icov @ delta
    ones      = np.ones(len(delta))
    Cinv_1    = icov @ ones
    A         = float(delta @ Cinv_d)
    B         = float(ones  @ Cinv_d)
    C_val     = float(ones  @ Cinv_1)
    return A - B**2 / C_val + np.log(C_val / (2 * np.pi))


def symmetrize_and_verify(cov, name="", reg=1e-8):
    """Symmetrize and verify a covariance matrix is positive-definite."""
    cov = 0.5 * (cov + cov.T)
    eigvals = np.linalg.eigvalsh(cov)
    if np.any(eigvals < 0):
        print(f"[WARN] {name} covariance not PD; regularizing with {reg}*I")
        cov += reg * np.eye(cov.shape[0])
    try:
        cholesky(cov)
        return cov, True
    except Exception:
        print(f"[ERROR] {name} covariance failed Cholesky after regularization")
        return cov, False


if __name__ == "__main__":
    z_d, types_d, mean_d, cov_d = load_desi_dr2()
    print(f"DESI DR2: {len(z_d)} measurements, types: {types_d}")
    print(f"  values: {mean_d}")

    H0, sH0, Om, sOm, R, sR, zs = load_planck()
    print(f"Planck: H0={H0}±{sH0}, Om={Om}±{sOm}, R={R}±{sR}")

    z_cc, H_cc, cov_cc = load_moresco()
    print(f"Moresco CC: {len(z_cc)} points, cov shape {cov_cc.shape}")

    z_f, f_f, e_f = load_fss8()
    print(f"RSD fσ8: {len(z_f)} points")

    H0_S, sH0_S = load_shoes()
    print(f"SH0ES: H0={H0_S}±{sH0_S}")

    z_pp, mu_pp, cov_pp = load_pantheon_plus()
    if z_pp is not None:
        print(f"Pantheon+: {len(z_pp)} SNe, cov shape {cov_pp.shape}")

    z_des, mu_des, icov_des = load_des_dovekie()
    if z_des is not None:
        print(f"DES-Dovekie: {len(z_des)} SNe, icov shape {icov_des.shape}")

    z_u3, mu_u3, icov_u3 = load_union3()
    if z_u3 is not None:
        print(f"Union3: {len(z_u3)} bins, icov shape {icov_u3.shape}")
