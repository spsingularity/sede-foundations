"""
SEDE joint MCMC pipeline using emcee.

Runs the full Bayesian fit with:
  - DESI DR2 BAO (13 measurements, full covariance)
  - Pantheon+ SN Ia (1580 SN, full STAT+SYS cov) — if data available
  - Planck CMB compressed (H0, Omega_m, R)
  - SH0ES H0 prior
  - Moresco CC (15 points, full Moresco-2020 cov)
  - fσ8 RSD (18 points)
  Cross-validation: DES-Dovekie, Union3

Parameters: (Omega_m, H0, gamma, rd, M_B)
"""

import numpy as np
import emcee
import os, time

from . import friedmann as fr
from . import data_loader as dl
from . import likelihood as lk
from .theory import w_DE_algebraic

# ─── Grid setup for pre-tabulation ──────────────────────────────────────────

OMEGA_M_GRID = np.linspace(0.20, 0.45, 26)
GAMMA_GRID   = np.concatenate([
    np.linspace(0.60, 3.0, 25),    # fine grid through theory value
    np.linspace(3.1, 10.0, 15),    # coarser at large gamma (→ΛCDM limit)
])
# sigma8 axis for the RUNNING-gamma SEDE-H model (replaces the gamma axis there:
# f_sat = Sigma_S(sigma8*D)/Sigma_S(sigma8) is sigma8-dependent, gamma_eff runs).
SIGMA8_GRID = np.linspace(0.60, 0.95, 12)
# Lean z-grid for the running table: dense where DE matters (z<12, the ODE
# region), sparse-log beyond (E=matter+rad analytic, only the integrated
# distance to z*~1090 needs coverage). Cuts n_z ~2x vs Z_GRID_TABLE.
Z_GRID_RUN = np.unique(np.concatenate([
    np.linspace(0.0, 5.0, 350),
    np.linspace(5.02, 12.0, 50),
    np.geomspace(12.1, 1500.0, 70),
]))
Z_GRID_TABLE = np.concatenate([
    np.linspace(0, 5, 500),
    np.geomspace(5.01, 1200, 500),
])
Z_GRID_TABLE = np.sort(np.unique(Z_GRID_TABLE))

from ._paths import CACHE_DIR   # layout-agnostic: <root>/cache in both dev-hub and release trees
CACHE_PATH       = os.path.join(CACHE_DIR, 'sede_table_v3.npz')
CACHE_PATH_LCDM  = os.path.join(CACHE_DIR, 'lcdm_table_v1.npz')
CACHE_PATH_SEDEH = os.path.join(CACHE_DIR, 'sedeH_table_v1.npz')
CACHE_PATH_SEDEH_RUN = os.path.join(CACHE_DIR, 'sedeH_run_table_v1.npz')
CACHE_PATH_COUSIN = CACHE_PATH_SEDEH   # naming-unification alias: 'cousin' == 'sedeH' (renamed model)

_CACHE_BY_MODEL = {'sede': CACHE_PATH, 'lcdm': CACHE_PATH_LCDM,
                   'sedeH': CACHE_PATH_SEDEH, 'cousin': CACHE_PATH_SEDEH}


def build_or_load_table(force_rebuild=False, model='sede'):
    """Build/load the (Omega_m x gamma x z) E(z), D_C(z) table for `model`."""
    cache = _CACHE_BY_MODEL[model]
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache) and not force_rebuild:
        print(f"Loading cached {model.upper()} table from {cache}")
        return fr.load_table(cache)
    print(f"Building {model.upper()} E(z), D_C(z) tables...")
    t0 = time.time()
    table = fr.build_table(OMEGA_M_GRID, GAMMA_GRID, Z_GRID_TABLE, cache, model=model)
    print(f"Table built in {time.time()-t0:.1f}s")
    return table


def build_or_load_table_running(force_rebuild=False, p=1.0):
    """
    Build/load the (Omega_m x sigma8 x z) table for the RUNNING-gamma SEDE-H
    model with entropy weight p (default p=1, the extensive/physical weight,
    Theorem 4B). The sigma8 axis is stored under 'gamma_grid' so the standard
    interpolators/likelihood work unchanged (the driver passes sigma8 in the
    structural-coordinate slot).
    """
    cache = CACHE_PATH_SEDEH_RUN
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache) and not force_rebuild:
        print(f"Loading cached SEDE-H-RUNNING table from {cache}")
        return fr.load_table(cache)
    print("Building SEDE-H-RUNNING table (exact structural entropy, p=%.3f)..." % p)
    from .gamma_computation import build_sigma_S_interp
    sigma_S = build_sigma_S_interp(p=p)
    t0 = time.time()
    table = fr.build_table_running(OMEGA_M_GRID, SIGMA8_GRID, sigma_S,
                                   z_grid=Z_GRID_RUN, cache_path=cache)
    print(f"Table built in {time.time()-t0:.1f}s")
    return table


def load_all_data():
    """Assemble all dataset dictionaries."""
    data = {}

    # DESI DR2
    z_d, types_d, mean_d, cov_d = dl.load_desi_dr2()
    cov_d, ok = dl.symmetrize_and_verify(cov_d, "DESI DR2")
    data['desi'] = (z_d, types_d, mean_d, cov_d)

    # Planck
    data['planck'] = dl.load_planck()

    # SH0ES
    data['shoes'] = dl.load_shoes()

    # Moresco
    z_cc, H_cc, cov_cc = dl.load_moresco()
    cov_cc, ok = dl.symmetrize_and_verify(cov_cc, "Moresco CC")
    data['moresco'] = (z_cc, H_cc, cov_cc)

    # fσ8
    data['fss8'] = dl.load_fss8()

    # SN Ia (optional — requires downloaded data)
    z_pp, mu_pp, cov_pp = dl.load_pantheon_plus()
    if z_pp is not None:
        cov_pp, ok = dl.symmetrize_and_verify(cov_pp, "Pantheon+")
        # Pre-compute Cholesky factor: saves O(n^3) at every MCMC step → O(n^2)
        try:
            L_pp = np.linalg.cholesky(cov_pp)
        except Exception:
            L_pp = None
        data['pantheon']   = (z_pp, mu_pp, cov_pp)
        data['pantheon_L'] = L_pp

    z_des, mu_des, icov_des = dl.load_des_dovekie()
    if z_des is not None:
        data['des'] = (z_des, mu_des, icov_des)

    z_u3, mu_u3, icov_u3 = dl.load_union3()
    if z_u3 is not None:
        data['union3'] = (z_u3, mu_u3, icov_u3)

    return data


# ─── LCDM baseline ───────────────────────────────────────────────────────────

def log_posterior_lcdm(theta_lcdm, data, interp_E_lcdm, interp_DC_lcdm):
    """
    LCDM log-posterior for comparison.
    theta_lcdm = (Omega_m, H0, rd, M_B).

    interp_E_lcdm / interp_DC_lcdm MUST come from a model='lcdm' table
    (E(z)=E_LCDM, gamma-independent). The gamma slot is then a dummy: we
    pass gamma=1.5 only so the (Omega_m, gamma, z) interpolation point is
    in-grid. (The old code embedded gamma=100, which extrapolated off the
    GAMMA_GRID and produced a ~10%-wrong E(z) → garbage LCDM chi^2.)
    """
    Omega_m, H0, rd, M_B, sigma8 = theta_lcdm
    theta_sede = (Omega_m, H0, 1.5, rd, M_B, sigma8)  # gamma is a no-op for the LCDM table
    return lk.log_posterior(theta_sede, data, interp_E_lcdm, interp_DC_lcdm)


# ─── MCMC runner ─────────────────────────────────────────────────────────────

def run_mcmc(data, interp_E, interp_DC,
             n_walkers=32, n_steps=1500, n_burnin=400,
             gamma_fixed=None,
             use_shoes=True, use_pantheon=True,
             label="SEDE"):
    """
    Run emcee MCMC.
    If gamma_fixed is not None, gamma is held fixed (reduces parameter space).

    Returns (sampler, flat_samples, param_names).
    """
    # Starting point. Params: (Omega_m, H0, gamma, rd, M_B, sigma8).
    # M_B ≈ -19.43 from Pantheon+ with H0=68.5, r_d=155; sigma8 starts at the
    # RSD-preferred ~0.78 (free parameter now, no longer pinned to Planck 0.811).
    theta0 = np.array([0.311, 67.4, 1.50, 147.0, -19.43, 0.78])
    pnames = ['Omega_m', 'H0', 'gamma', 'rd', 'M_B', 'sigma8']
    scale_of = {'Omega_m': 0.01, 'H0': 0.5, 'gamma': 0.05,
                'rd': 0.5, 'M_B': 0.02, 'sigma8': 0.02}

    if gamma_fixed is not None:
        # Fix gamma — sample (Omega_m, H0, rd, M_B, sigma8)
        theta0 = np.delete(theta0, 2)
        pnames.pop(2)

        def log_prob(theta):
            Om, H0, rd, MB, s8 = theta
            th_full = (Om, H0, gamma_fixed, rd, MB, s8)
            return lk.log_posterior(th_full, data, interp_E, interp_DC,
                                     use_shoes=use_shoes,
                                     use_pantheon=use_pantheon)
    else:
        def log_prob(theta):
            return lk.log_posterior(theta, data, interp_E, interp_DC,
                                     use_shoes=use_shoes,
                                     use_pantheon=use_pantheon)

    ndim = len(theta0)
    scales = np.array([scale_of[n] for n in pnames])
    rng   = np.random.default_rng(42)
    p0    = theta0 + scales * rng.standard_normal((n_walkers, ndim))

    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_prob)

    print(f"\n[{label}] Burn-in ({n_burnin} steps, {n_walkers} walkers)...")
    t0 = time.time()
    state = sampler.run_mcmc(p0, n_burnin, progress=True)
    sampler.reset()

    print(f"[{label}] Production ({n_steps} steps)...")
    sampler.run_mcmc(state, n_steps, progress=True)
    dt = time.time() - t0

    acc = np.mean(sampler.acceptance_fraction)
    print(f"[{label}] Done in {dt:.0f}s. Acceptance fraction: {acc:.3f}")
    if not 0.15 < acc < 0.75:
        print(f"  WARNING: acceptance fraction {acc:.3f} outside [0.15, 0.75]")

    flat = sampler.get_chain(flat=True)
    return sampler, flat, pnames


# ─── Results summary ─────────────────────────────────────────────────────────

def summarize_results(flat_samples, pnames, data, interp_E, interp_DC,
                      gamma_fixed=None):
    """Print parameter constraints and chi^2 breakdown."""
    medians = np.median(flat_samples, axis=0)
    lo, hi  = np.percentile(flat_samples, [16, 84], axis=0)

    print("\n" + "=" * 55)
    print("SEDE MCMC RESULTS")
    print("=" * 55)
    for i, (name, med, l, h) in enumerate(zip(pnames, medians, lo, hi)):
        print(f"  {name:10s} = {med:.4f}  +{h-med:.4f}  -{med-l:.4f}")

    Om_best = medians[pnames.index('Omega_m')]
    H0_best = medians[pnames.index('H0')]
    rd_best = medians[pnames.index('rd')]
    MB_best = medians[pnames.index('M_B')]
    s8_best = medians[pnames.index('sigma8')]
    gam_best = gamma_fixed if gamma_fixed is not None else medians[pnames.index('gamma')]

    theta_best = (Om_best, H0_best, gam_best, rd_best, MB_best, s8_best)

    # Algebraic EOS prediction (Theorem 5)
    w0_pred = w_DE_algebraic(Om_best)
    print(f"\n  w_DE(0) algebraic (Theorem 5) = {w0_pred:.4f}")
    print(f"  DESI DR2 w_0                  = -0.838 ± 0.055")
    print(f"  Tension                        = {abs(w0_pred - (-0.838))/0.055:.2f} sigma")

    # chi^2 breakdown
    print("\nChi^2 breakdown:")
    sectors = [
        ('DESI DR2',  lambda: lk.chi2_desi(theta_best, data, interp_E, interp_DC),  13),
        ('Planck CMB', lambda: lk.chi2_planck(theta_best, data, interp_E, interp_DC), 3),
        ('SH0ES',      lambda: lk.chi2_shoes(theta_best, data),                       1),
        ('Moresco CC', lambda: lk.chi2_moresco(theta_best, data, interp_E, interp_DC), 15),
        ('fσ8 RSD',    lambda: lk.chi2_fss8(theta_best, data, interp_E, interp_DC),   18),
    ]
    if data.get('pantheon') and data['pantheon'][0] is not None:
        sectors.append(('Pantheon+', lambda: lk.chi2_pantheon(theta_best, data, interp_E, interp_DC), 1580))

    total_chi2 = 0.0
    total_dof  = 0
    for name, fn, dof in sectors:
        try:
            c2 = fn()
            total_chi2 += c2
            total_dof  += dof
            print(f"  {name:12s}: chi^2={c2:.1f}  dof={dof}  chi^2/dof={c2/dof:.2f}")
        except Exception as e:
            print(f"  {name:12s}: ERROR {e}")

    n_params = len(pnames)
    print(f"\n  Total chi^2 = {total_chi2:.1f}  dof={total_dof-n_params}")
    print(f"  AIC = chi^2 + 2k = {total_chi2 + 2*n_params:.1f}")
    print(f"  BIC = chi^2 + k*ln(n) = {total_chi2 + n_params*np.log(total_dof):.1f}")


if __name__ == "__main__":
    print("Building table...")
    table = build_or_load_table()
    interp_E, interp_DC = fr.make_interpolators(table)
    fr.verify_table(table)

    print("\nLoading data...")
    data = load_all_data()

    print("\nRunning SEDE MCMC (gamma fixed to theory)...")
    from .gamma_computation import compute_gamma
    gamma_th, _, _ = compute_gamma()
    print(f"gamma_theory = {gamma_th:.4f}")

    sampler, flat, pnames = run_mcmc(
        data, interp_E, interp_DC,
        n_walkers=32, n_steps=1500, n_burnin=400,
        gamma_fixed=gamma_th,
        label="SEDE",
    )
    summarize_results(flat, pnames, data, interp_E, interp_DC,
                      gamma_fixed=gamma_th)
