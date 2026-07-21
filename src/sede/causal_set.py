"""
Causal-set test of SEDE's dof-COUNTING premise (the counting half of the §8.5 postulate).

SEDE's Δ is fixed not by the horizon's quantum STATE (chaos/SYK settles that — V60) but by how
the number of horizon degrees of freedom scales with radius: area N ∝ R^{d-2} (⟹ Δ=0,
Bekenstein–Hawking) vs spatial volume N ∝ R^{d-1} (⟹ Δ=1, SEDE). Causal-set quantum gravity —
Lorentz-invariantly discrete spacetime, points Poisson-sprinkled at Planck density — is the natural
setting for this question, and it contains BOTH scalings:

  BULK    — # elements in the spatial ball of radius R          ~ R^{d-1}  (spatial volume)
  HORIZON — # causal LINKS (covering relations) across radius R ~ R^{d-2}  (area; the causal-set
            area law that recovers Bekenstein–Hawking, Sorkin and collaborators)

Worked in 2+1D (d=3): BULK ~ R², HORIZON ~ R¹; the one-power gap is the Δ=1-vs-Δ=0 distinction.
FINDING (honest): the *canonical* horizon-entropy object (links) is AREA-law; SEDE's Δ=1 requires
identifying s_grav with the BULK count instead — a non-default choice even in volume-friendly QG.
"""
import numpy as np


def sprinkle(N, L, T, rng):
    """Poisson sprinkle N points into 2+1D Minkowski box t∈[-T,T], (x,y)∈[-L,L]²."""
    return np.column_stack([rng.uniform(-T, T, N),
                            rng.uniform(-L, L, N),
                            rng.uniform(-L, L, N)])


def sprinkle_dS(N, L, eta0, U, rng):
    """Poisson sprinkle N points into a 2+1D de Sitter patch in CONFORMAL coordinates (η,x,y),
    η∈[-U,-eta0] (η<0), with the dS proper-volume measure dV ∝ a³ dη d²x, a=-1/(Hη) ⟹ density
    ∝ |η|^{-3}. de Sitter is conformally flat, so the causal relation is IDENTICAL to Minkowski's
    in these coordinates (Δη>0, Δη²>Δx²+Δy²) — the whole point: the causal-set area law transfers
    to dS unchanged. Returns columns [η, x, y]; sample η by inverse-CDF of |η|^{-3}."""
    v = rng.uniform(0, 1, N)
    # inverse CDF of p(u)∝u^{-3} on [eta0,U]:  u = (eta0^{-2} - v(eta0^{-2}-U^{-2}))^{-1/2}
    u = (eta0**-2 - v * (eta0**-2 - U**-2))**-0.5
    eta = -u
    x = rng.uniform(-L, L, N)
    y = rng.uniform(-L, L, N)
    return np.column_stack([eta, x, y])


def bulk_count(P, R, t0=0.0, dt=0.4):
    """# elements in the spatial disk r<R within a thin time slab |t-t0|<dt (∝ spatial volume)."""
    m = np.abs(P[:, 0] - t0) < dt
    r = np.hypot(P[m, 1], P[m, 2])
    return int(np.sum(r < R))


def horizon_links(P, R, tau_cut, shell):
    """# causal LINKS (covering relations) straddling the radius-R circle (r_in<R<r_out).
    Links are short (proper time < tau_cut), so both endpoints and any intermediate element lie
    within `shell` of r=R; restrict to shell points and vectorise the emptiness test (∝ area)."""
    r = np.hypot(P[:, 1], P[:, 2])
    Pn = P[np.abs(r - R) < shell]
    rn = r[np.abs(r - R) < shell]
    order = np.argsort(Pn[:, 0])
    Pn, rn = Pn[order], rn[order]
    t, x, y = Pn[:, 0], Pn[:, 1], Pn[:, 2]
    inside = rn < R
    n = len(Pn)
    nlink = 0
    for i in range(n):
        js = np.arange(i + 1, n)
        if len(js) == 0:
            break
        dtj = t[js] - t[i]
        win = dtj <= tau_cut
        js, dtj = js[win], dtj[win]
        if len(js) == 0:
            continue
        sep2 = (x[js] - x[i])**2 + (y[js] - y[i])**2
        cj = js[(dtj * dtj > sep2) & (inside[js] != inside[i])]   # related & straddling
        for j in cj:
            zt = (t > t[i]) & (t < t[j])
            if not zt.any():
                nlink += 1
                continue
            zi = np.where(zt)[0]
            iz = (t[zi] - t[i])**2 > (x[zi] - x[i])**2 + (y[zi] - y[i])**2
            zj = (t[j] - t[zi])**2 > (x[j] - x[zi])**2 + (y[j] - y[zi])**2
            if not np.any(iz & zj):                               # empty interval ⟹ link
                nlink += 1
    return nlink


def power_law(counts, Rs):
    """Power-law exponent of counts vs Rs (least-squares in log-log)."""
    c, R = np.asarray(counts, float), np.asarray(Rs, float)
    g = c > 0
    return float(np.polyfit(np.log(R[g]), np.log(c[g]), 1)[0])


def run_experiment(rho=6.0, L=10.0, T=10.0, Rs=(3, 4, 5, 6, 7), n_real=6, seed=0):
    """Sprinkle and measure bulk (∝volume) and horizon-link (∝area) scalings. Returns dict."""
    rng = np.random.default_rng(seed)
    N = int(rho * (2 * L)**2 * (2 * T))
    tau_cut = 1.6 * rho**(-1 / 3)
    shell = 2.0 * tau_cut
    Rs = list(Rs)
    bulk = np.zeros(len(Rs))
    horiz = np.zeros(len(Rs))
    for _ in range(n_real):
        P = sprinkle(N, L, T, rng)
        for k, R in enumerate(Rs):
            bulk[k] += bulk_count(P, R)
            horiz[k] += horizon_links(P, R, tau_cut, shell)
    bulk /= n_real
    horiz /= n_real
    e_bulk, e_horiz = power_law(bulk, Rs), power_law(horiz, Rs)
    return dict(N=N, rho=rho, Rs=Rs, bulk=bulk.tolist(), horizon=horiz.tolist(),
                exp_bulk=e_bulk, exp_horizon=e_horiz, gap=e_bulk - e_horiz)


def run_experiment_dS(rho=200.0, L=10.0, eta0=1.0, U=8.0, t0=-2.0,
                      Rs=(3, 4, 5, 6, 7), n_real=6, seed=0):
    """Same area-vs-volume measurement in de SITTER (conformal coords, dS measure). Because dS is
    conformally flat the causal relation is unchanged, so the link area-law and bulk volume-law
    BOTH transfer to dS — the dS cosmological horizon does NOT prefer volume. Returns same dict."""
    rng = np.random.default_rng(seed)
    # mean local density near the analysis band sets tau_cut/shell
    rho_local = rho / abs(t0)**3
    N = int(rho * (2 * L)**2 * 0.5 * (eta0**-2 - U**-2))   # = rho * spatial-area * ∫|η|^{-3}dη
    tau_cut = 1.6 * rho_local**(-1 / 3)
    shell = 2.0 * tau_cut
    Rs = list(Rs)
    bulk = np.zeros(len(Rs))
    horiz = np.zeros(len(Rs))
    for _ in range(n_real):
        P = sprinkle_dS(N, L, eta0, U, rng)
        for k, R in enumerate(Rs):
            bulk[k] += bulk_count(P, R, t0=t0, dt=0.5)
            horiz[k] += horizon_links(P[np.abs(P[:, 0] - t0) < 2.0], R, tau_cut, shell)
        del P
    bulk /= n_real
    horiz /= n_real
    e_bulk, e_horiz = power_law(bulk, Rs), power_law(horiz, Rs)
    return dict(geometry="de Sitter", N=N, rho=rho, Rs=Rs, bulk=bulk.tolist(),
                horizon=horiz.tolist(), exp_bulk=e_bulk, exp_horizon=e_horiz, gap=e_bulk - e_horiz)
