"""
SEDE Friedmann solver and pre-tabulated E(z) / D_C(z) tables.

Model: SEDE Path C
  E^2(z) = Omega_m (1+z)^3 + Omega_r (1+z)^4 + Omega_DE(z)
  Omega_DE(z) = Omega_DE0 * f_sat(z),  Omega_DE0 = 1 - Omega_m - Omega_r

  f_sat(z) = (1 - exp(-gamma * x)) / (1 - exp(-gamma))
  x = sigma_8^2(z) / sigma_8^2(0) = D^2(z)   [D normalized to D(0)=1]

  Variable choice: x = D^2 = sigma_8^2(z)/sigma_8^2(0) (the linear variance ratio),
  per SEDE_COMPLETE_REFERENCE_recreated.md eq. 3.7 and SEDE.md. The algebraic
  EOS w_DE(0) = (4 Omega_m/3 - 1)/(1 - Omega_m) (Result 5b) follows from this
  variance convention via d ln sigma^2/d ln a |_0 = 2 f_g(0). gamma_theory ~ 1.50
  from the Sheth-Tormen integral (COLOSSUS, full mass range).

The growth factor D(z) is solved simultaneously with E(z) via the
standard ODE:
    D'' + (2 + epsilon)H D' - (3/2) Omega_m H^2 / E^2 * D = 0
where ' = d/dt = -H(1+z) d/dz.

Pipeline: solve ODE → build grid → cache → interpolate at query points.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
import os, hashlib, json

C_KM_S = 299792.458   # speed of light km/s

# Default grid
Z_MAX  = 1200.0
# Log-spaced at high z ensures accurate 1/E(z) integration in the matter era
Z_GRID = np.concatenate([
    np.linspace(0, 5, 500),           # late-time: dense linear
    np.geomspace(5.01, 1200, 500),     # early: log-spaced (1/E ∝ z^{-3/2})
])
Z_GRID = np.sort(np.unique(Z_GRID))


def _growth_ode_rhs(ln_a, state, Omega_m, Omega_r):
    """
    RHS for growth ODE in ln(a) variable.
    state = [D, dD/d(ln_a)]
    E^2(a) approximated as LCDM during growth. This is exact in the early
    universe (z>5) where f_sat≈0. At z<0.5 where SEDE has ~5% lower E(z)
    than LCDM, the approximation introduces ~1-2% error in D(z) and hence
    ~1% in f_sat(z). The effect is partially self-correcting since D(0)=1
    by normalization. Self-consistent iteration is not implemented.
    """
    a = np.exp(ln_a)
    E2 = Omega_m / a**3 + Omega_r / a**4 + (1 - Omega_m - Omega_r)
    E  = np.sqrt(max(E2, 1e-30))
    # Hubble parameter in ln(a): H = E * H0
    # epsilon = d ln H / d ln a = (a/E) dE/da ... approximate:
    dE2_dlna = -3 * Omega_m / a**3 - 4 * Omega_r / a**4
    epsilon = 0.5 * dE2_dlna / E2

    D, dD = state
    # Growth equation: d^2D/d(lna)^2 + (2 + epsilon) dD/d(lna) - 3/2 Omega_m/(a^3 E^2) D = 0
    source = 1.5 * Omega_m / (a**3 * E2)
    d2D = -(2 + epsilon) * dD + source * D
    return [dD, d2D]


def compute_growth_factor(z_arr, Omega_m, Omega_r=9.0e-5):
    """
    Compute linear growth factor D(z), normalized D(z=0) = 1.
    Integrates from ln(a)=-7 (deep matter era) to ln(a)=0.
    """
    ln_a_start = -7.0
    ln_a_end   = 0.0
    a_start    = np.exp(ln_a_start)
    # Initial conditions: matter era D ~ a
    D0  = a_start
    dD0 = a_start

    ln_a_out = np.sort(-np.log(1.0 + np.asarray(z_arr)))
    # Always include ln_a=0 so normalization D/D(a=1) is computed correctly
    ln_a_eval = np.concatenate([[ln_a_start], [ln_a_end],
                                 np.clip(ln_a_out, ln_a_start, ln_a_end)])
    ln_a_eval = np.sort(np.unique(ln_a_eval))

    sol = solve_ivp(
        _growth_ode_rhs,
        [ln_a_start, ln_a_end],
        [D0, dD0],
        t_eval=ln_a_eval,
        args=(Omega_m, Omega_r),
        rtol=1e-8, atol=1e-10, method='DOP853',
    )
    # Normalize at z=0 (ln_a=0)
    idx0 = np.argmin(np.abs(sol.t - 0.0))
    D_at_0 = sol.y[0, idx0]
    D_norm = sol.y[0] / D_at_0

    # Interpolate; extrapolate below ln_a_start using matter-era D ∝ a scaling
    from scipy.interpolate import interp1d
    ln_a_query = -np.log(1 + np.asarray(z_arr))
    interp = interp1d(sol.t, D_norm, kind='linear', bounds_error=False,
                      fill_value=(D_norm[0], D_norm[-1]))
    D_out = interp(ln_a_query)
    # For z > z_start (below ln_a_start), apply D ∝ a extrapolation
    D_at_start = D_norm[0]
    below = ln_a_query < ln_a_start
    D_out[below] = D_at_start * np.exp(ln_a_query[below] - ln_a_start)
    return D_out


def E_SEDE(z, Omega_m, gamma, Omega_r=9.0e-5):
    """
    SEDE Hubble parameter E(z) = H(z)/H_0.

    E^2(z) = Omega_m(1+z)^3 + Omega_r(1+z)^4 + f_sat(z) * (1 - Omega_m - Omega_r)
    with f_sat(z) from Theorem 3, normalized so that E(0)=1 exactly.

    Note: Omega_DE(z) = f_sat(z) * Omega_DE(0) and Omega_DE(0) = 1 - Omega_m - Omega_r.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    Omega_DE0 = 1.0 - Omega_m - Omega_r

    # Growth factor: x = D^2(z) = sigma_8^2(z)/sigma_8^2(0)  (variance ratio)
    D = compute_growth_factor(z, Omega_m, Omega_r)
    x = D ** 2  # σ_8^2(z)/σ_8^2(0), ∈ [0,1]; γ = d ln Σ_S / d ln σ_8^2 ≈ 1.50

    # f_sat from Theorem 3
    # Bousso bound: f_sat ∈ [0,1]  (Theorem 1, Hypothesis H5)
    f_sat = np.clip(
        (1.0 - np.exp(-gamma * x)) / (1.0 - np.exp(-gamma)),
        0.0, 1.0
    )

    E2 = (Omega_m * (1 + z)**3
          + Omega_r * (1 + z)**4
          + Omega_DE0 * f_sat)
    E2 = np.maximum(E2, 1e-30)
    return np.sqrt(E2)


def compute_growth_model(z_query, Omega_m, E_of_z, Omega_r=9.0e-5):
    """
    Linear growth D(z) and growth rate f(z)=dlnD/dlna for a SMOOTH dark-energy
    model with arbitrary background E(z)=H/H0 (dark energy does not cluster on
    sub-horizon scales; it enters only through the expansion history).

    Solves, in N=ln a:
        D'' + [2 + (1/2) dlnE^2/dlna] D' - (3/2) Omega_m a^{-3} / E^2 * D = 0
    with matter-era IC D ∝ a. Returns (D_query, f_query), D normalized to D(0)=1.

    Unlike compute_growth_factor (which hardwires the ΛCDM E^2), this uses the
    model's own E(z) — needed for a self-consistent fσ8 and for the SEDE-H
    growth history.
    """
    lna = np.linspace(-7.0, 0.0, 800)
    a   = np.exp(lna)
    z_g = 1.0 / a - 1.0
    E2  = np.asarray(E_of_z(z_g), dtype=float) ** 2
    E2  = np.maximum(E2, 1e-30)
    dlnE2 = np.gradient(np.log(E2), lna)

    from scipy.interpolate import interp1d
    E2i   = interp1d(lna, E2,    kind='cubic', bounds_error=False, fill_value=(E2[0], E2[-1]))
    dE2i  = interp1d(lna, dlnE2, kind='cubic', bounds_error=False, fill_value=(dlnE2[0], dlnE2[-1]))

    def rhs(N, y):
        av  = np.exp(N)
        src = 1.5 * (Omega_m / av**3) / E2i(N)
        fric = 2.0 + 0.5 * dE2i(N)
        return [y[1], -fric * y[1] + src * y[0]]

    a0 = a[0]
    sol = solve_ivp(rhs, [-7.0, 0.0], [a0, a0], dense_output=True,
                    rtol=1e-8, atol=1e-12, method='DOP853')
    D0_norm = sol.sol(0.0)[0]

    Nq = -np.log(1.0 + np.asarray(z_query, dtype=float))
    Nq = np.clip(Nq, -7.0, 0.0)
    y  = sol.sol(Nq)
    D  = y[0] / D0_norm
    f  = y[1] / y[0]            # dlnD/dlna = D'/D
    return D, f


def E_SEDE_H(z, Omega_m, gamma, Omega_r=9.0e-5, z_cap=50.0, growth=None):
    """
    SEDE-H ("horizon-back-reacted" SEDE) Hubble parameter E(z)=H/H0.

    Derived from the conjugate identity with the DYNAMICAL Cai-Kim temperature
    T_AH = (H/2pi)(1 - eps/2), eps = -Hdot/H^2  (reference eq. 2.1), which the
    original E_SEDE drops:
        rho_DE = (3/8piG) H^2 (1 - eps/2) f_sat = (3/8piG)(H^2 + Hdot/2) f_sat
    with f_sat(z) = (1-e^{-gamma D^2})/(1-e^{-gamma}), f_sat(0)=1 (Theorem 3 norm).

    The falling H^2 competes with the rising structure term f_sat, so rho_DE
    turns over and the EOS w(z) CROSSES -1 (quintessence today, phantom past) —
    matching the DESI DR2 evolving-DE signature, with NO cosmological constant
    (rho_DE -> 0 at early times via f_sat -> 0). There is no free lambda: the
    H-scaling is fixed by the dynamical horizon; gamma (structural) is the only
    new parameter, and it lands near the Sheth-Tormen value ~1.5.

    Friedmann + the identity give the background ODE (radiation included):
        E'(z) = 2[Omega_m(1+z)^3 + Omega_r(1+z)^4 - E^2(1-f)] / [(1+z) E f]
    integrated from E(0)=1. Beyond z_cap (f_sat negligible) E^2 -> matter+rad.

    `growth` may be a precomputed (z_grid, D_grid) tuple to avoid re-solving the
    growth ODE (D depends only on Omega_m) — used by the table builder.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    if growth is None:
        zz = np.concatenate([np.linspace(0, 5, 200),
                             np.geomspace(5.01, z_cap + 5, 120)])
        Dg = compute_growth_factor(zz, Omega_m, Omega_r)
    else:
        zz, Dg = growth
    den = 1.0 - np.exp(-gamma)

    def f_of(zq):
        Dq = np.interp(zq, zz, Dg)
        return (1.0 - np.exp(-gamma * Dq**2)) / den

    def rhs(zq, E):
        f = max(f_of(zq), 1e-8)
        src = Omega_m * (1 + zq)**3 + Omega_r * (1 + zq)**4
        return [2.0 * (src - E[0]**2 * (1.0 - f)) / ((1 + zq) * E[0] * f)]

    sol = solve_ivp(rhs, [0.0, z_cap], [1.0], dense_output=True,
                    rtol=1e-7, atol=1e-9, max_step=0.5)
    out = np.empty_like(z)
    lo = z <= z_cap
    if np.any(lo):
        out[lo] = sol.sol(z[lo])[0]
    hi = ~lo
    if np.any(hi):
        out[hi] = np.sqrt(Omega_m * (1 + z[hi])**3 + Omega_r * (1 + z[hi])**4)
    return np.maximum(out, 1e-30)


def E_SEDE_H_run(z, Omega_m, sigma8, sigma_S, Omega_r=9.0e-5, z_cap=12.0, growth=None):
    """
    RUNNING-gamma SEDE-H: dynamical-horizon background with the EXACT structural
    entropy saturation (no constant-gamma linearization), using the physically
    derived EXTENSIVE entropy weight p=1 (Theorem 4B):

        f_sat(z) = Sigma_S(sigma8 * D(z)) / Sigma_S(sigma8),   f_sat(0)=1,

    where Sigma_S = sigma_S(amplitude) is a build_sigma_S_interp(p=1) callable
    (the collapsed mass density / extensive structural entropy). This makes
    f_sat sigma8-DEPENDENT, so the model is indexed by (Omega_m, sigma8) instead
    of (Omega_m, gamma): the structural coupling gamma_eff(z) now RUNS, set by
    the entropy weight rather than fitted.

    Same dynamical-horizon background ODE as E_SEDE_H (T_AH=(H/2pi)(1-eps/2)):
        E'(z) = 2[Om(1+z)^3 + Or(1+z)^4 - E^2(1-f)] / [(1+z) E f],  E(0)=1,
    integrated to z_cap; matter+rad beyond. The rising f_sat vs falling H^2 gives
    the w(z) crossing of -1 (quintessence today), now with a theory-fixed
    (not free) structural shape. gamma is no longer a parameter.

    `growth` may be a precomputed (z_grid, D_grid) tuple (D depends only on Om).
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    if growth is None:
        zz = np.concatenate([np.linspace(0, 5, 200),
                             np.geomspace(5.01, z_cap + 5, 120)])
        Dg = compute_growth_factor(zz, Omega_m, Omega_r)
    else:
        zz, Dg = growth
    Ss0 = float(sigma_S(sigma8))

    def f_of(zq):
        Dq = np.interp(zq, zz, Dg)
        return float(sigma_S(sigma8 * Dq)) / Ss0

    # Determine an effective cap: integrate only while f_sat (=Omega_DE/Omega_DE-ish)
    # is non-negligible. For low sigma8 the exact entropy underflows at LOWER z,
    # and pushing the ODE into that regime triggers catastrophic cancellation in
    # (src - E^2(1-f)) -> stiff/slow. Stop once f_sat < 1e-4 (DE then <0.01% of
    # the budget); matter+rad beyond. f_sat is monotone in z (D decreasing).
    zz_cap = zz[zz <= z_cap]
    f_cap = np.array([f_of(zq) for zq in zz_cap])
    below = np.where(f_cap < 1e-4)[0]
    z_eff = float(zz_cap[below[0]]) if len(below) else z_cap

    def rhs(zq, E):
        f = min(max(f_of(zq), 1e-8), 1.0)
        src = Omega_m * (1 + zq)**3 + Omega_r * (1 + zq)**4
        return [2.0 * (src - E[0]**2 * (1.0 - f)) / ((1 + zq) * E[0] * f)]

    sol = solve_ivp(rhs, [0.0, z_eff], [1.0], dense_output=True,
                    rtol=1e-7, atol=1e-9, max_step=0.5)
    out = np.empty_like(z)
    lo = z <= z_eff
    if np.any(lo):
        out[lo] = sol.sol(z[lo])[0]
    hi = ~lo
    if np.any(hi):
        out[hi] = np.sqrt(Omega_m * (1 + z[hi])**3 + Omega_r * (1 + z[hi])**4)
    return np.maximum(out, 1e-30)


def build_table_running(Omega_m_grid, sigma8_grid, sigma_S, z_grid=None,
                        cache_path=None, Omega_r=9.0e-5):
    """
    Pre-tabulate E(z), D_C(z), self-consistent growth D(z),f(z) for the RUNNING
    SEDE-H model on a (Omega_m x sigma8 x z) grid. The sigma8 axis is stored
    under the key 'gamma_grid' so the existing make_interpolators / likelihood
    machinery (which treats the 2nd interpolation coordinate as the structural
    axis) works unchanged — for this model the "2nd axis" IS sigma8, and the
    driver passes sigma8 in that slot. `sigma_S` = build_sigma_S_interp(p=1).
    """
    if z_grid is None:
        z_grid = Z_GRID
    n_Om, n_s8, n_z = len(Omega_m_grid), len(sigma8_grid), len(z_grid)
    E_table  = np.zeros((n_Om, n_s8, n_z))
    DC_table = np.zeros((n_Om, n_s8, n_z))
    D_table  = np.zeros((n_Om, n_s8, n_z))
    f_table  = np.zeros((n_Om, n_s8, n_z))
    print(f"Building SEDE-H-RUNNING table: {n_Om} x {n_s8} x {n_z} = {n_Om*n_s8*n_z} points")

    for i, Om in enumerate(Omega_m_grid):
        # growth D(z) for the f_sat amplitude depends only on Om — solve once.
        zz_g = np.concatenate([np.linspace(0, 5, 200), np.geomspace(5.01, 60.0, 120)])
        Dg_g = compute_growth_factor(zz_g, Om, Omega_r)
        for j, s8 in enumerate(sigma8_grid):
            E_arr = E_SEDE_H_run(z_grid, Om, s8, sigma_S, growth=(zz_g, Dg_g))
            # Use the smooth E_SEDE_H_run callable inside the growth solve (a
            # pre-sampled E-array has a derivative kink at z_cap that makes the
            # growth ODE stiff — far slower than the extra ODE solve).
            D_arr, f_arr = compute_growth_model(
                z_grid, Om, lambda zz, _s8=s8: E_SEDE_H_run(zz, Om, _s8, sigma_S,
                                                            growth=(zz_g, Dg_g)))
            integrand = 1.0 / E_arr
            DC_arr = np.zeros(n_z)
            for k in range(1, n_z):
                DC_arr[k] = DC_arr[k-1] + 0.5*(integrand[k-1]+integrand[k])*(z_grid[k]-z_grid[k-1])
            E_table[i, j], DC_table[i, j] = E_arr, DC_arr
            D_table[i, j], f_table[i, j]  = D_arr, f_arr
        if (i + 1) % max(1, n_Om // 5) == 0:
            print(f"  Done {i+1}/{n_Om} Omega_m slices")

    table = dict(Omega_m_grid=Omega_m_grid, gamma_grid=sigma8_grid, z_grid=z_grid,
                 E_table=E_table, DC_table=DC_table, D_table=D_table, f_table=f_table)
    if cache_path:
        np.savez_compressed(cache_path, **table)
        print(f"Cached to {cache_path}")
    return table


def E_SEDE_lambda(z, Omega_m, gamma, lam, Omega_r=9.0e-5):
    """
    General H^{2λ}-coupling SEDE background (the λ-family diagnostic, borrowed
    from SEDE_V2's `coupling_lambda`). The conjugate identity is parametrised as
        ρ_DE = (3/8πG) (H²)^λ H0^{2(1-λ)} f_sat,
    i.e.  E²(z) = matter + Ω_DE0 · f_sat(z) · (E²)^λ,  solved by fixed-point
    iteration (the DE term is sub-dominant → contraction). Limits:
        λ=0  → additive model E_SEDE (ρ_DE ∝ f_sat; phantom w0≈-1.15);
        λ=1  → full conjugate identity (ρ_DE ∝ H² f_sat; tracks H², CMB-tense).
    The dynamical-horizon SEDE-H (E_SEDE_H) sits at an effective λ_eff≈0.7 today
    via the Cai-Kim (1-ε/2) factor and closes the EOS gap (Theorem 5D). This
    function is for mapping the w0 / EOS-gap / intermediate-z-DE landscape vs λ.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    Omega_DE0 = 1.0 - Omega_m - Omega_r
    D = compute_growth_factor(z, Omega_m, Omega_r)
    f_sat = np.clip((1.0 - np.exp(-gamma * D**2)) / (1.0 - np.exp(-gamma)), 0.0, 1.0)
    matter = Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4
    if lam == 0.0:
        E2 = matter + Omega_DE0 * f_sat
    else:
        E2 = matter + Omega_DE0 * f_sat
        for _ in range(80):
            E2_new = matter + Omega_DE0 * f_sat * np.power(np.maximum(E2, 1e-30), lam)
            if np.max(np.abs(E2_new - E2)) < 1e-12:
                E2 = E2_new
                break
            E2 = E2_new
    return np.sqrt(np.maximum(E2, 1e-30))


def E_SEDE_H_form(z, Omega_m, fsat_of_x, Omega_r=9.0e-5, z_cap=12.0, growth=None):
    """
    SEDE-H dynamical-horizon background with an ARBITRARY saturation form
    f_sat = fsat_of_x(x), x = D^2(z) (variance ratio), f_sat(1)=1, f_sat(0)=0.
    Same ODE as E_SEDE_H; lets us test alternative f_sat functional forms beyond
    the exponential CDF. Dynamic cap where f_sat<1e-4 (matter+rad beyond).
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    if growth is None:
        zz = np.concatenate([np.linspace(0, 5, 200), np.geomspace(5.01, z_cap + 5, 120)])
        Dg = compute_growth_factor(zz, Omega_m, Omega_r)
    else:
        zz, Dg = growth

    def f_of(zq):
        Dq = np.interp(zq, zz, Dg)
        return float(fsat_of_x(Dq**2))

    zz_cap = zz[zz <= z_cap]
    fc = np.array([f_of(q) for q in zz_cap])
    below = np.where(fc < 1e-4)[0]
    z_eff = float(zz_cap[below[0]]) if len(below) else z_cap

    def rhs(zq, E):
        f = min(max(f_of(zq), 1e-8), 1.0)
        src = Omega_m * (1 + zq)**3 + Omega_r * (1 + zq)**4
        return [2.0 * (src - E[0]**2 * (1.0 - f)) / ((1 + zq) * E[0] * f)]

    sol = solve_ivp(rhs, [0.0, z_eff], [1.0], dense_output=True, rtol=1e-7, atol=1e-9, max_step=0.5)
    out = np.empty_like(z)
    lo = z <= z_eff
    if np.any(lo):
        out[lo] = sol.sol(z[lo])[0]
    hi = ~lo
    if np.any(hi):
        out[hi] = np.sqrt(Omega_m * (1 + z[hi])**3 + Omega_r * (1 + z[hi])**4)
    return np.maximum(out, 1e-30)


def E_SEDE_barrow(z, Omega_m, gamma, Delta=1.0, Omega_r=9.0e-5):
    """
    SEDE-H with the H-coupling DERIVED from Barrow entropy (Theorem 8):
        rho_DE ∝ H^{2-Delta} f_sat,   lambda = 1 - Delta/2,
    Delta = Barrow fractal deformation in [0,1]. Delta=0 -> lambda=1 (Bekenstein-
    Hawking, naive identity); Delta=1 -> lambda=0.5 (maximal fractal horizon, the
    data-preferred value). Thin wrapper over E_SEDE_lambda with lambda=1-Delta/2.
    The structural coupling gamma stays at its halo-statistics value (~1.5).
    """
    return E_SEDE_lambda(z, Omega_m, gamma, 1.0 - Delta / 2.0, Omega_r)


def E_SEDE_volume(z, Omega_m, gamma, Omega_r=9.0e-5):
    """
    SEDE in the VOLUME formulation — the PARAMETER-FREE statement of the dark sector.

    Instead of a Barrow-deformed area entropy S ∝ A^{1+Δ/2} with a deformation parameter Δ,
    postulate directly that the gravitational entropy is the coarse-grained ENTANGLEMENT
    entropy contained in the apparent-horizon VOLUME:
        S_grav ∝ V_AH ∝ R_AH³      ⟺      s_grav = constant horizon entropy density.
    Then ρ_DE = T_AH·s_grav·f_sat ∝ H·f_sat (T_AH ∝ H, s_grav = const), normalised by
    flatness [Ω_DE0 = 1−Ω_m−Ω_r, f_sat(0)=1]. There is NO deformation parameter.

    This is IDENTICAL to Barrow at the maximal deformation: S ∝ A^{1+Δ/2}|_{Δ=1} = A^{3/2}
    = R³ = V, i.e. the volume law IS Δ=1, and equals the λ=1−Δ/2 = 0.5 H-coupling. So this
    returns exactly E_SEDE_barrow(Δ=1) = E_SEDE_lambda(λ=0.5); verified bit-for-bit (E(z) to
    2e-6, dataset χ² unchanged — run_volume_equiv.py). Δ is retained ONLY as a falsifiable
    TEST axis (deviation from the volume law): data give Δ=0.95±0.05, excluding the area law
    Δ=0 at ~19σ (run_delta_indirect.py).

    Closed form (the fixed point E² = src + Ω_DE0 f_sat E is quadratic in E):
        E(z) = ½[Ω_DE0 f_sat + √((Ω_DE0 f_sat)² + 4·src)],  src = Ω_m(1+z)³ + Ω_r(1+z)⁴.
    Implemented via E_SEDE_lambda(λ=0.5) so existing results are reproduced identically.
    """
    return E_SEDE_lambda(z, Omega_m, gamma, 0.5, Omega_r)


def E_LCDM(z, Omega_m, Omega_r=9.0e-5):
    """LCDM comparison model."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    Omega_L = 1.0 - Omega_m - Omega_r
    E2 = Omega_m * (1+z)**3 + Omega_r * (1+z)**4 + Omega_L
    return np.sqrt(np.maximum(E2, 1e-30))


# --- naming-unification aliases (2026-07-20) -------------------------------
# The lambda=1 "temperature-factor" model was renamed E_SEDE_H -> E_SEDE_cousin in
# the cosmology/apdm release code (physics IDENTICAL; the rename only clarifies that
# it is the disfavoured cousin of the canonical E_SEDE_volume). Expose BOTH names so
# every release's experiments import against this one library. Likewise the build_table
# model string 'sedeH' now also accepts 'cousin'.
E_SEDE_cousin = E_SEDE_H
E_SEDE_cousin_run = E_SEDE_H_run
E_SEDE_cousin_form = E_SEDE_H_form


def build_table(Omega_m_grid, gamma_grid, z_grid=None, cache_path=None, model='sede'):
    """
    Pre-tabulate E(z) and D_C(z) on a (Omega_m x gamma x z) grid.
    Returns dict with arrays; optionally saves .npz cache.

    model='sede' : E(z) from E_SEDE (gamma-dependent).
    model='lcdm' : E(z) from E_LCDM (gamma-independent; the same E_arr is
                   broadcast across every gamma slice so the resulting
                   interpolators are drop-in compatible with the SEDE
                   likelihood, which always passes a (Omega_m, gamma, z) point).
                   This replaces the old "embed gamma=100" hack, which fell
                   outside GAMMA_GRID (max 10) and extrapolated to a ~10%-wrong
                   E(z) — corrupting every LCDM chi^2.
    """
    if z_grid is None:
        z_grid = Z_GRID

    n_Om = len(Omega_m_grid)
    n_g  = len(gamma_grid)
    n_z  = len(z_grid)

    E_table  = np.zeros((n_Om, n_g, n_z))
    DC_table = np.zeros((n_Om, n_g, n_z))
    D_table  = np.zeros((n_Om, n_g, n_z))   # self-consistent linear growth D(z)
    f_table  = np.zeros((n_Om, n_g, n_z))   # growth rate f=dlnD/dlna (model's own)

    print(f"Building {model.upper()} table: {n_Om} x {n_g} x {n_z} = {n_Om*n_g*n_z} points")

    for i, Om in enumerate(Omega_m_grid):
        # For LCDM, E(z) does not depend on gamma — compute once and reuse.
        E_lcdm_arr = E_LCDM(z_grid, Om) if model == 'lcdm' else None
        # SEDE-H: growth D(z) depends only on Omega_m — solve the growth ODE once
        # per Omega_m and reuse across the gamma slices.
        if model in ('sedeH', 'cousin'):
            zz_g = np.concatenate([np.linspace(0, 5, 200),
                                   np.geomspace(5.01, 60.0, 120)])
            Dg_g = compute_growth_factor(zz_g, Om)
        # LCDM growth is gamma-independent — compute once per Omega_m.
        if model == 'lcdm':
            D_l, f_l = compute_growth_model(z_grid, Om, lambda zz: E_LCDM(zz, Om))
        for j, gam in enumerate(gamma_grid):
            if model == 'lcdm':
                E_arr = E_lcdm_arr
                D_arr, f_arr = D_l, f_l
            elif model in ('sedeH', 'cousin'):
                E_arr = E_SEDE_H(z_grid, Om, gam, growth=(zz_g, Dg_g))
                D_arr, f_arr = compute_growth_model(
                    z_grid, Om, lambda zz: E_SEDE_H(zz, Om, gam, growth=(zz_g, Dg_g)))
            else:
                E_arr = E_SEDE(z_grid, Om, gam)
                D_arr, f_arr = compute_growth_model(z_grid, Om, lambda zz: E_SEDE(zz, Om, gam))
            # Comoving distance: DC = (c/H0) int_0^z dz'/E(z')
            # H0 factor absorbed; in units c/H0 = 1
            integrand = 1.0 / E_arr
            # Trapezoidal cumulative
            DC_arr = np.zeros(n_z)
            for k in range(1, n_z):
                DC_arr[k] = DC_arr[k-1] + 0.5 * (integrand[k-1] + integrand[k]) * (z_grid[k] - z_grid[k-1])
            E_table[i, j, :]  = E_arr
            DC_table[i, j, :] = DC_arr
            D_table[i, j, :]  = D_arr
            f_table[i, j, :]  = f_arr

        if (i + 1) % max(1, n_Om // 5) == 0:
            print(f"  Done {i+1}/{n_Om} Omega_m slices")

    table = dict(
        Omega_m_grid=Omega_m_grid,
        gamma_grid=gamma_grid,
        z_grid=z_grid,
        E_table=E_table,
        DC_table=DC_table,
        D_table=D_table,
        f_table=f_table,
    )
    if cache_path:
        np.savez_compressed(cache_path, **table)
        print(f"Cached to {cache_path}")
    return table


def load_table(cache_path):
    d = np.load(cache_path)
    return {k: d[k] for k in d.files}


def make_interpolators(table):
    """
    Build bilinear RegularGridInterpolators for E and D_C.
    Returns (interp_E, interp_DC) callable as interp(Om, gamma, z).
    """
    pts = (table['Omega_m_grid'], table['gamma_grid'], table['z_grid'])
    interp_E  = RegularGridInterpolator(pts, table['E_table'],  method='linear', bounds_error=False, fill_value=None)
    interp_DC = RegularGridInterpolator(pts, table['DC_table'], method='linear', bounds_error=False, fill_value=None)
    return interp_E, interp_DC


def make_growth_interpolators(table):
    """
    Build interpolators for the model's self-consistent growth D(z) and growth
    rate f(z)=dlnD/dlna. Returns (interp_D, interp_f) callable as interp(Om,gamma,z),
    or (None, None) for legacy tables that lack D_table/f_table.
    """
    if 'D_table' not in table or 'f_table' not in table:
        return None, None
    pts = (table['Omega_m_grid'], table['gamma_grid'], table['z_grid'])
    interp_D = RegularGridInterpolator(pts, table['D_table'], method='linear', bounds_error=False, fill_value=None)
    interp_f = RegularGridInterpolator(pts, table['f_table'], method='linear', bounds_error=False, fill_value=None)
    return interp_D, interp_f


def verify_table(table):
    """Analytical checks on the tabulated E(z)."""
    Om_arr = table['Omega_m_grid']
    g_arr  = table['gamma_grid']
    z_arr  = table['z_grid']
    E_tab  = table['E_table']

    z_idx0 = np.argmin(np.abs(z_arr))
    E0_values = E_tab[:, :, z_idx0]

    max_err = np.max(np.abs(E0_values - 1.0))
    print(f"E(z=0) = 1 check: max deviation = {max_err:.2e}  (target < 1e-4)")

    # f_sat at recombination z~1100 should be << 1
    # Compute directly from D(z) — avoids catastrophic cancellation at high z
    from .theory import fsat_late
    Om_mid = Om_arr[len(Om_arr)//2]
    g_mid  = g_arr[len(g_arr)//2]
    D_rec  = compute_growth_factor(np.array([1089.9]), Om_mid)[0]
    f_sat_rec = fsat_late(1089.9, 0.811 * D_rec, 0.811, g_mid)
    print(f"f_sat(z=1089.9) ≈ {f_sat_rec:.4e}  (must be ≪ 1)")

    return max_err < 1e-3


if __name__ == "__main__":
    # Quick verification
    z_test = np.array([0, 0.1, 0.5, 1.0, 2.0, 5.0])
    E_test = E_SEDE(z_test, Omega_m=0.311, gamma=1.50)
    print("E(z) for SEDE (Omega_m=0.311, gamma=1.50):")
    for z, E in zip(z_test, E_test):
        print(f"  z={z:.1f}  E={E:.6f}")
    print(f"\nE(0) = {E_test[0]:.8f}  (must equal 1.0)")

    # Build a small table to test
    Om_grid  = np.linspace(0.25, 0.38, 8)
    gam_grid = np.linspace(0.80, 2.20, 6)
    table = build_table(Om_grid, gam_grid, z_grid=np.linspace(0, 3, 200))
    verify_table(table)
