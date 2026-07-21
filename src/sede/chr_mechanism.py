"""
CHR — Critical Horizon Response (Theorem 14).

Unifies the three candidate mechanisms for SEDE's structure-sourcing into ONE
relaxational order-parameter dynamics for the horizon entropy state:

  (A) criticality / susceptibility  — why a tiny driver gates an O(1) activation
  (C) gravitational backreaction     — the structure variance IS the control field
  (E) interacting dark energy        — the energy-ledger reading of the relaxation

The existing logistic f_sat (Theorems 3, 12) is the ADIABATIC EQUILIBRIUM of this
dynamics. CHR adds the response sector — susceptibility, fluctuation enhancement,
environment dependence — that distinguishes SEDE from any kinematic w(z), plus one
falsifiable identity:

  CHR IDENTITY (one number, three hats; cf. the Theorem-11 seam):
        ε_dep   =   1/χ   =   𝒮''   ≈ 1.5×10⁻⁷
     (deposited-entropy fraction) = (inverse susceptibility) = (proximity to spinodal)

i.e. the tiny structure-deposited entropy and the large near-critical amplification
are the SAME number — the horizon sits within ε_dep of the area↔volume spinodal,
driven there by the generalised second law (self-organised criticality, not tuning).

See run_chr_experiments.py for the six predictions and THEOREM_14 in theory.py.
"""
import numpy as np
from .friedmann import compute_growth_factor, compute_growth_model, E_SEDE_lambda

GAMMA = 1.4964          # γ = horizon-deposited entropy weight (Theorem 4C), p = 5/3
LAMBDA = 0.5            # λ = 1 − Δ/2, Δ = 1 (Theorems 8, 9)


# ----------------------------------------------------------------------------
# Order parameter, control field, susceptibility (the Landau triplet)
# ----------------------------------------------------------------------------
def f_eq(x, gamma=GAMMA):
    """Equilibrium order parameter φ_eq (= the existing f_sat). x = D² = σ²-ratio."""
    return np.clip((1.0 - np.exp(-gamma * x)) / (1.0 - np.exp(-gamma)), 0.0, 1.0)


def susceptibility(x, gamma=GAMMA):
    """
    χ_x(x) ≡ dφ_eq/dx — the DERIVED response of the order parameter to the structure-variance
    driver x = D². This is the slope of the fixed kernel, an O(1) number (≈1.93 at x=0, ≈0.43 at
    x=1 for γ=1.5) — no free parameter. It is the susceptibility that enters the background w(z)
    lock (P2). Distinct from χ_J = ∂φ/∂J = χ_x·(dx/dJ) ∼ 1/ε_dep ∼ 10⁷, the response to the
    microscopic deposited-entropy FIELD, which appears only in the OPTIONAL near-critical layer
    (P3–P5) — see chi_J() and Theorem 14.
    """
    return gamma * np.exp(-gamma * x) / (1.0 - np.exp(-gamma))


def chi_J(Omega_m=0.30, **kw):
    """
    χ_J ≡ ∂φ/∂J ∼ 1/ε_dep — the SECOND susceptibility, response to the deposited-entropy field J.
    It is merely the ratio of the O(1) variance to the 10⁻⁷ entropy fraction (χ_J = χ_x·dx/dJ),
    NOT a fitted parameter and NOT needed for the background w(z) result; it sets the amplitude of
    the optional critical-response predictions only. Returns ~7×10⁶.
    """
    return 1.0 / deposited_entropy_fraction(Omega_m, **kw)


def control_variance(z, Omega_m=0.30, Omega_r=9e-5):
    """
    Control field of the transition = structure variance ratio x ≡ σ²(z)/σ²(0) = D²(z).
    This is the (normalised) kinematical backreaction Q_𝒟 ∝ ⟨θ²⟩−⟨θ⟩² of forming
    structure (mechanism C); the volume-law bound (Theorem 9) supplies the Buchert
    closure. f_sat = f(x) is therefore backreaction-driven by construction.
    """
    D = compute_growth_factor(np.atleast_1d(z), Omega_m, Omega_r)
    return D**2


# ----------------------------------------------------------------------------
# (A) The CHR identity:  ε_dep = 1/χ = 𝒮''  (one number, three hats)
# ----------------------------------------------------------------------------
def deposited_entropy_fraction(Omega_m=0.30, f_coll=0.5, sigma_v_kms=300.0):
    """
    ε_dep = ΔS_bind / S_AH — the fraction of the horizon (Gibbons–Hawking) entropy
    represented by the structure-released binding entropy deposited at T_AH (Theorem 4C).

    Order of magnitude (Gibbons–Hawking T·S = E for the area-law reservoir, so the
    entropy fraction equals the binding-energy density fraction):
        ε_dep ≈ Ω_m · f_coll · ⟨σ_v²⟩/c²
    with collapsed fraction f_coll ~ 0.5 and virial velocities σ_v ~ 300 km/s. This is
    the referee's "binding entropy ≪ horizon entropy" (~10⁻⁷) — and the SAME 6-orders
    short as the binding-ENERGY-vs-dark-energy gap, by construction.
    """
    c = 299792.458
    return Omega_m * f_coll * (sigma_v_kms / c) ** 2


def chr_identity(Omega_m=0.30, **kw):
    """
    Returns the CHR identity numbers:
      ε_dep   — deposited-entropy fraction (tiny field J)
      chi     — susceptibility required for φ = O(1):  χ ≈ 1/ε_dep
      Spp     — entropy-landscape curvature 𝒮'' ≈ ε_dep (proximity to the spinodal)
    The point: ε_dep = 1/χ = 𝒮'' is ONE number. The horizon is self-organised to within
    ε_dep of the area↔volume spinodal by the GSL (Δ→1 attractor, Theorem 9).

    HONEST STATUS (checkpoint review): χ ≡ 1/ε_dep and 𝒮'' ≡ ε_dep are DEFINITIONAL here, so
    the "CHR identity" ε_dep·χ = 1 holds by construction — it is NOT an independent numerical
    result. The physical content is the POSTULATE that the horizon sits within ε_dep of the
    spinodal (near-criticality S''≈ε_dep, Theorem 14, flagged OPEN). Whether that near-critical
    proximity is generic (SOC attractor) or tuned is the open question; it is assumed here, not
    derived. Downstream "the drive is amplified to O(1)" is therefore conditional on this postulate.
    """
    eps_dep = deposited_entropy_fraction(Omega_m, **kw)
    chi = 1.0 / eps_dep          # definitional (near-critical susceptibility ≡ 1/ε_dep)
    Spp = eps_dep                # definitional (spinodal proximity ≡ ε_dep) — the Thm-14 postulate
    return dict(eps_dep=eps_dep, chi=chi, Spp=Spp)


# ----------------------------------------------------------------------------
# Unified EOS:  1+w = (1/3)[ 2λε − (χ/φ)·2x·g ]   (A·χ, C·x, growth-lock g)
# ----------------------------------------------------------------------------
def eos_decomposition(z, Omega_m=0.30, gamma=GAMMA, lam=LAMBDA, Omega_r=9e-5):
    """
    The CHR reading of Theorem 13: adiabatic elimination of the order parameter gives
        1 + w_DE(z) = (1/3) [ 2λ ε(z)  −  (χ/φ)·2x·g ],
    where ε = −dlnH/dlna (expansion), and the phantom term factorises into
        χ   = susceptibility            (mechanism A)
        x   = structure variance        (mechanism C, backreaction driver)
        g   = dlnD/dlna = growth rate   (the RSD/WL observable that LOCKS w to growth)

    Returns z, and dicts with: w_direct (from ρ_DE), w_chr (the factorised formula),
    and the individual terms. The phantom term uses ONLY growth observables (x, g) —
    that is the w(z)↔σ²(z) lock no kinematic w(z) obeys.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    E = E_SEDE_lambda(z, Omega_m, gamma, lam, Omega_r)
    eps = (1.0 + z) * np.gradient(np.log(E), z)                       # ε = −dlnH/dlna

    # growth: D(z) and g = dlnD/dlna for the SEDE background itself
    D, g = compute_growth_model(z, Omega_m, lambda zz: E_SEDE_lambda(np.atleast_1d(zz),
                                Omega_m, gamma, lam, Omega_r), Omega_r)
    x = D**2
    phi = np.clip(f_eq(x, gamma), 1e-30, 1.0)
    chi = susceptibility(x, gamma)

    expansion = (2.0 * lam * eps) / 3.0                              # pushes w > −1
    # entropy-production term = (1/3) dln f/dlna = (1/3)(χ/φ)·dx/dlna, dx/dlna = 2 x g
    entropy = (1.0 / 3.0) * (chi / phi) * (2.0 * x * g)              # pushes w < −1
    w_chr = -1.0 + expansion - entropy

    rho = np.maximum(E**2 - Omega_m * (1 + z)**3 - Omega_r * (1 + z)**4, 1e-12)
    w_direct = -1.0 + (1.0 / 3.0) * (1 + z) * np.gradient(np.log(rho), z)

    return z, dict(w_direct=w_direct, w_chr=w_chr, expansion=expansion,
                   entropy=entropy, eps=eps, x=x, g=g, chi=chi, phi=phi, E=E)


def transition_redshift(Omega_m=0.30, gamma=GAMMA, Omega_r=9e-5, zmax=4.0):
    """z* where the order parameter crosses ½ (the critical/transition epoch; ≈1.2)."""
    z = np.linspace(0.0, zmax, 2000)
    phi = f_eq(control_variance(z, Omega_m, Omega_r), gamma)
    return float(z[np.argmin(np.abs(phi - 0.5))])


# ----------------------------------------------------------------------------
# (E) Interacting-DE energy exchange  Q_int = ρ_DE φ̇/φ
# ----------------------------------------------------------------------------
def interaction_rate(z, Omega_m=0.30, gamma=GAMMA, lam=LAMBDA, Omega_r=9e-5):
    """
    Q_int / (H ρ_DE) = (1/φ) dφ/d(lna) = (χ/φ)·2x·g — the dimensionless dark-sector
    coupling, FIXED by the structure-growth rate (not a free β). Coupled quintessence
    in PPF-legible form (Theorem 6B / V43 handle the w=−1 crossing).
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    D, g = compute_growth_model(z, Omega_m, lambda zz: E_SEDE_lambda(np.atleast_1d(zz),
                                Omega_m, gamma, lam, Omega_r), Omega_r)
    x = D**2
    phi = np.clip(f_eq(x, gamma), 1e-30, 1.0)
    return (susceptibility(x, gamma) / phi) * (2.0 * x * g)
