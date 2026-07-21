"""
Compute gamma from the Sheth-Tormen mass function via COLOSSUS.
Theorem 4 of the SEDE framework.

gamma = +d ln Sigma_S / d ln sigma_8^2 = +(1/2) d ln Sigma_S / d ln sigma_8
where Sigma_S = integral_{M_min}^{M_max} M^{5/3} dn/dM dM

Sign: Sigma_S INCREASES with sigma_8 (more halos at higher clustering),
so d ln Sigma_S / d ln sigma_8 > 0, giving gamma > 0. This positive
gamma enters f_sat as exp(-gamma*x), a DECREASING exponential with
x = D^2(z) = sigma_8^2(z)/sigma_8^2(0).

Convention: x = D^2 = sigma_8^2/sigma_8^2(0) (variance ratio), per
SEDE_COMPLETE_REFERENCE_recreated.md eq. 3.7. The coupling is
gamma = d ln Sigma_S / d ln sigma_8^2 = (1/2) d ln Sigma_S / d ln sigma_8 ~ 1.50
(COLOSSUS, full mass range 10^10-10^16 M_sun, Sheth-Tormen).
"""

import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import mass_function

# Mass range and reference cosmology
LOG_M_MIN = 10.0   # log10(M / M_sun)
LOG_M_MAX = 16.0
N_MASS_BINS = 500
MASS_ARRAY = np.logspace(LOG_M_MIN, LOG_M_MAX, N_MASS_BINS)  # M_sun/h

# Fiducial cosmology (Planck 2018)
FIDUCIAL_PARAMS = dict(
    H0=67.4, Om0=0.315, Ob0=0.049,
    sigma8=0.811, ns=0.965, Tcmb0=2.725,
)


def _set_cosmo(sigma8_val):
    """Set COLOSSUS cosmology with modified sigma8."""
    params = dict(FIDUCIAL_PARAMS)
    params['sigma8'] = sigma8_val
    name = f"sede_s8_{sigma8_val:.4f}"
    try:
        cosmology.addCosmology(name, params)
    except Exception:
        pass
    cosmology.setCosmology(name)


def compute_entropy_weighted_integral(sigma8_val):
    """
    Compute Sigma_S(sigma8) = int M^{5/3} dn/dM dM.
    Uses Sheth-Tormen mass function from COLOSSUS.
    """
    _set_cosmo(sigma8_val)
    mf = mass_function.massFunction(
        MASS_ARRAY, z=0.0,
        mdef='fof', model='sheth99',
        q_out='dndlnM',
    )
    # dn/dM from dn/dlnM: dn/dM = (dn/dlnM) / M
    dndM = mf / MASS_ARRAY
    integrand = MASS_ARRAY ** (5.0 / 3.0) * dndM
    # Trapezoidal integration in log-M
    sigma_S = np.trapezoid(integrand, MASS_ARRAY)
    return sigma_S


def compute_gamma(delta_ln_s8=0.02):
    """
    Compute gamma = +d ln Sigma_S / d ln sigma_8 via finite differences.

    Uses two-point centered difference:
        gamma = +[ln Sigma_S(s8+) - ln Sigma_S(s8-)] / [2 delta ln sigma_8]

    Convention: f_sat(z) = (1-exp(-gamma*D^2(z)))/(1-exp(-gamma)) with x=D^2.
    Therefore gamma = d ln Sigma_S / d ln sigma_8^2 = (1/2) d ln Sigma_S / d ln sigma_8
    (derivative w.r.t. variance). At the fiducial Sheth-Tormen calibration:
    gamma ≈ 1.50 (Theorem 4).
    """
    s8_fid = FIDUCIAL_PARAMS['sigma8']
    s8_hi = s8_fid * np.exp(+delta_ln_s8)
    s8_lo = s8_fid * np.exp(-delta_ln_s8)

    Ss_hi  = compute_entropy_weighted_integral(s8_hi)
    Ss_fid = compute_entropy_weighted_integral(s8_fid)
    Ss_lo  = compute_entropy_weighted_integral(s8_lo)

    # gamma = d ln Sigma_S / d ln sigma_8^2 = (1/2) d ln Sigma_S / d ln sigma_8
    # x = D^2(z) = sigma_8^2(z)/sigma_8^2(0), so gamma is the slope w.r.t. ln(sigma_8^2)
    dlnSs_dlns8 = (np.log(Ss_hi) - np.log(Ss_lo)) / (2 * delta_ln_s8)
    gamma = +0.5 * dlnSs_dlns8  # ~ 1.50

    _set_cosmo(s8_fid)
    return gamma, dlnSs_dlns8, Ss_fid


def compute_gamma_at_sigma8(s8, delta_ln_s8=0.02):
    """
    Local structural slope gamma = (1/2) d ln Sigma_S / d ln sigma8 evaluated at
    an ARBITRARY clustering amplitude s8 (= sigma8 at some epoch). gamma is the
    f_sat exponent at that pivot amplitude.
    """
    s8_hi = s8 * np.exp(+delta_ln_s8)
    s8_lo = s8 * np.exp(-delta_ln_s8)
    Ss_hi = compute_entropy_weighted_integral(s8_hi)
    Ss_lo = compute_entropy_weighted_integral(s8_lo)
    return 0.5 * (np.log(Ss_hi) - np.log(Ss_lo)) / (2 * delta_ln_s8)


def gamma_under_growth(D_func, z_arr=(0.0, 0.5, 1.0, 2.0), s8_0=0.811):
    """
    Re-derive gamma_theory using a SPECIFIED linear growth history D(z)
    (e.g. SEDE-H's self-consistent growth) instead of the default ΛCDM one.

    The structural coupling is gamma = d ln Sigma_S / d ln sigma8^2 evaluated at
    the pivot amplitude. Since sigma8(z) = sigma8_0 * D(z) and D(0)=1 by
    normalization, the z=0 value gamma_theory is INVARIANT under the growth
    choice (it is an amplitude-space derivative of the z=0 halo abundance) — this
    is reference Tier-2 Proof 5. The growth model only changes the *running*
    gamma_eff(z) at z>0 (where sigma8 differs between models). The f_sat formula
    uses the z=0 pivot, so gamma_theory is growth-independent.

    Returns dict: {'gamma0', 'z', 'gamma_eff', 'sigma8_z'}.
    """
    z_arr = np.asarray(z_arr, dtype=float)
    D = np.asarray(D_func(z_arr), dtype=float)
    s8_z = s8_0 * D
    gamma_eff = np.array([compute_gamma_at_sigma8(s) for s in s8_z])
    gamma0 = compute_gamma_at_sigma8(s8_0)
    return {'gamma0': gamma0, 'z': z_arr, 'gamma_eff': gamma_eff, 'sigma8_z': s8_z}


def build_sigma_S_interp(p=5.0/3.0, s8_min=0.04, s8_max=0.95, n=28):
    """
    Build an interpolator for the entropy-weighted halo abundance
        Sigma_S^{(p)}(sigma8) = ∫ M^p (dn/dM)(M; sigma8) dM
    as a function of the clustering amplitude sigma8, for an arbitrary entropy
    weight exponent p (p=5/3 = gravitational self-energy; p=1 = mass/extensive).

    Returns a callable Sigma_S(sigma8). Used by the RUNNING-gamma f_sat:
        f_sat(z) = Sigma_S(sigma8 * D(z)) / Sigma_S(sigma8),   D(0)=1,
    which builds the redshift-running of gamma_eff in exactly (no constant-gamma
    linearization). The log-slope (1/2) dlnSigma_S/dlnsigma8 IS gamma_eff at that
    amplitude.
    """
    from scipy.interpolate import interp1d
    s8_grid = np.linspace(s8_min, s8_max, n)
    lnSs = np.empty(n)
    for i, s in enumerate(s8_grid):
        _set_cosmo(s)
        mf = mass_function.massFunction(MASS_ARRAY, z=0.0, mdef='fof',
                                        model='sheth99', q_out='dndlnM')
        dndM = mf / MASS_ARRAY
        lnSs[i] = np.log(np.trapezoid(MASS_ARRAY**p * dndM, MASS_ARRAY))
    _set_cosmo(FIDUCIAL_PARAMS['sigma8'])
    f = interp1d(np.log(s8_grid), lnSs, kind='cubic', fill_value='extrapolate')
    return lambda s: np.exp(f(np.log(np.clip(s, 1e-3, None))))


def fsat_structural(D, sigma8, sigma_S):
    """
    RUNNING-gamma saturation function from the exact structural entropy:
        f_sat = Sigma_S(sigma8 * D) / Sigma_S(sigma8),  clipped to [0,1].
    `D` is the linear growth factor (D(0)=1), `sigma_S` a build_sigma_S_interp().

    KEY RESULT (2026-06-23): with the gravitational-self-energy weight p=5/3,
    gamma_eff RUNS from ~1.5 at z=0 to ~20 at high z (massive halos are
    exponentially rare), so f_sat collapses far too steeply and the joint fit is
    REJECTED (DESI+SN chi^2 ~ 1461 vs LCDM ~1400; w(1)~-1.7). The data instead
    prefer a MILDER entropy weight p~1 (entropy ∝ mass), giving gamma_eff(z=0)~0.3
    and a fit that BEATS LCDM (chi^2 ~ 1396.5). So the constant-gamma_data~1.0 of
    the joint MCMC is consistent with extensive (p~1) entropy, NOT with the steep
    p=5/3 running — that is the resolution of the gamma_data-vs-theory gap.
    """
    return np.clip(np.asarray(sigma_S(sigma8 * D)) / sigma_S(sigma8), 0.0, 1.0)


def entropy_weight_scan(s8=0.811, delta_ln_s8=0.02, p_list=(2.0/3.0, 1.0, 4.0/3.0, 5.0/3.0)):
    """
    Scan the z=0 structural coupling gamma_eff(z=0) = (1/2) d ln Sigma_S^{(p)}/d ln sigma8
    over entropy-weight exponents p (Theorem 4B). Returns list of (p, gamma_eff).

    Key discriminants:
      p=2/3 (holographic area S∝A∝R^2):  gamma_eff < 0  -> EXCLUDED (f_sat runs
            the wrong way; DE would grow toward the past).
      p=1   (extensive entropy S∝N∝M):   gamma_eff ≈ +0.27 -> the physical weight
            (Theorem 4B); Sigma_S = collapsed mass density.
      p=5/3 (gravitational self-ENERGY): gamma_eff ≈ +1.50 -> the old Theorem-4
            value; rejected by the running fit (gamma_eff blows up to ~20 at high z).

    This makes explicit that the joint-fit preference gamma_data≈1 selects the
    extensive (entropy) weight, not the energy weight — see Theorem 4B.
    """
    out = []
    for p in p_list:
        Ss = build_sigma_S_interp(p=p)
        g = 0.5 * (np.log(Ss(s8 * np.exp(+delta_ln_s8)))
                   - np.log(Ss(s8 * np.exp(-delta_ln_s8)))) / (2 * delta_ln_s8)
        out.append((float(p), float(g)))
    return out


def gamma_eff_running(D_func, p=1.0, z_arr=(0.0, 0.5, 1.0, 2.0), s8_0=0.811,
                      delta_ln_s8=0.02):
    """
    Redshift-running structural coupling gamma_eff(z) = (1/2) d ln Sigma_S^{(p)}/d ln sigma8
    evaluated at the running amplitude sigma8(z) = s8_0 * D(z), for entropy weight p.

    With the physical extensive weight p=1 (Theorem 4B) this runs mildly,
    gamma_eff: ~0.27 (z=0) -> ~0.9 (z=2), bracketing the constant gamma_data≈1.0
    of the free-gamma joint MCMC. With p=5/3 it runs steeply (->~20), which the
    data reject. Returns dict {'z','sigma8_z','gamma_eff'}.
    """
    Ss = build_sigma_S_interp(p=p)
    z_arr = np.asarray(z_arr, dtype=float)
    s8_z = s8_0 * np.asarray(D_func(z_arr), dtype=float)
    g = np.array([
        0.5 * (np.log(Ss(s * np.exp(+delta_ln_s8))) - np.log(Ss(s * np.exp(-delta_ln_s8))))
        / (2 * delta_ln_s8) for s in s8_z])
    return {'z': z_arr, 'sigma8_z': s8_z, 'gamma_eff': g}


def collapsed_fraction(s8, p=1.0):
    """
    Sigma_S^{(p)}(sigma8). For p=1 this is the collapsed MASS density
    rho_coll(>M_min) = ∫ M (dn/dM) dM (Theorem 4B): the extensive-entropy
    structural budget. f_sat = collapsed_fraction(s8*D)/collapsed_fraction(s8).
    """
    Ss = build_sigma_S_interp(p=p)
    return float(Ss(s8))


def gamma_convergence_check():
    """Verify gamma is stable across different finite-difference step sizes."""
    steps = [0.01, 0.02, 0.05, 0.10]
    gammas = []
    print("Convergence of gamma with finite-difference step size:")
    print(f"  {'delta_ln_s8':>15}  {'gamma':>8}")
    for d in steps:
        g, _, _ = compute_gamma(delta_ln_s8=d)
        gammas.append(g)
        print(f"  {d:>15.3f}  {g:>8.4f}")
    return gammas


def gamma_mass_range_sensitivity():
    """Check sensitivity to mass range (cluster-only vs full)."""
    ranges = [
        (10, 16, "Full range [10^10, 10^16]"),
        (13, 16, "Cluster-only [10^13, 10^16]"),
        (10, 14, "Galaxy-scale [10^10, 10^14]"),
    ]
    s8_fid = FIDUCIAL_PARAMS['sigma8']
    results = []
    print("\nGamma sensitivity to mass range:")
    print(f"  {'Range':>32}  {'gamma':>8}")
    for lmin, lmax, label in ranges:
        global MASS_ARRAY
        MASS_ARRAY = np.logspace(lmin, lmax, N_MASS_BINS)
        g, _, _ = compute_gamma()
        results.append((label, g))
        print(f"  {label:>32}  {g:>8.4f}")
    # Restore
    MASS_ARRAY = np.logspace(LOG_M_MIN, LOG_M_MAX, N_MASS_BINS)
    return results


if __name__ == "__main__":
    print("Computing gamma from Sheth-Tormen mass function (COLOSSUS)...")
    gamma, dlnSs, Ss_fid = compute_gamma()
    print(f"\ngamma_theory = {gamma:.4f}")
    print(f"Expected:      ~1.50  (d ln Sigma_S / d ln sigma_8^2, x=D^2 convention)")
    print(f"d ln Sigma_S / d ln sigma_8 = {dlnSs:.4f}")
    print(f"Sigma_S (fiducial) = {Ss_fid:.4e} M_sun^{{2/3}} Mpc^{{-3}}")

    print("\n--- Convergence check ---")
    gamma_convergence_check()

    print("\n--- Mass range sensitivity ---")
    gamma_mass_range_sensitivity()

    # Final value for use in pipeline
    gamma_final, _, _ = compute_gamma(delta_ln_s8=0.02)
    print(f"\nAdopted gamma_theory = {gamma_final:.4f}")
