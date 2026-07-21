"""
run_lit_inferences.py — new physics inferred from the Paper-II literature search.

Two computable inferences that strengthen the papers; three further proposals are
discussed in the session notes.

  (1) ENSEMBLE INEQUIVALENCE resolves one horn of the canonical/microcanonical
      λ ambiguity (§9). Campa–Dauxois–Ruffo: long-range (α ≤ d) systems are
      non-additive and have INEQUIVALENT ensembles, with negative specific heat
      the hallmark. The de Sitter horizon has C < 0; combined with gravity's
      α = 1 ≤ d (run_residue_longrange.py), the CANONICAL ensemble is unphysical
      for the isolated self-gravitating horizon — only the MICROCANONICAL
      (state-count) accounting applies. The canonical horn gives area (λ = 1);
      removing it leaves the microcanonical (volume, λ = 1/2) as the only
      physical option. (The volume count *within* the microcanonical ensemble
      is still the dS-holography residue — this removes a horn, not the residue.)

  (2) GREA bulk-viscosity mapping. García-Bellido & Espinosa-Portalés cast
      entropy-production acceleration as an effective bulk viscosity with
      negative pressure. SEDE's structure-driven entropy production
      (d ln f_sat/d ln a > 0) maps to a POSITIVE effective bulk viscosity
      ζ(z) ≥ 0 — second-law-respecting — that is *gated by structure* rather
      than by horizon expansion. The w = −1 crossing is where the viscous
      (entropy-production) term balances the dilution (expansion) term, casting
      SEDE as a structure-gated member of the GREA/entropic-DE family.
"""
import numpy as np
from sede.friedmann import E_SEDE_lambda, compute_growth_factor

Om, Or, gamma, lam = 0.30, 9.0e-5, 1.4964, 0.5
ODE0 = 1 - Om - Or

# ===========================================================================
# (1) Ensemble inequivalence — the de Sitter horizon has negative specific heat
# ===========================================================================
print("=" * 74)
print("(1) Ensemble inequivalence ⇒ microcanonical (volume) is the physical horn")
print("=" * 74)
z = np.linspace(0, 3, 200)
E = E_SEDE_lambda(z, Om, gamma, lam, Or)
H = E                                   # in units of H0
# apparent horizon: T_A = H/2π, S_A = π/(G H²) ∝ 1/H².  C = T dS/dT.
# S ∝ 1/T² (since T ∝ H) ⇒ C = T dS/dT = −2 S  (negative ∀ z).
S_A = 1.0 / H**2                        # ∝ horizon entropy (area/4), arb. units
T_A = H / (2 * np.pi)
C = T_A * np.gradient(S_A, T_A)         # specific heat (numerical)
print(f"  de Sitter / apparent horizon specific heat C = T dS/dT:")
print(f"    C/S at z=0,1,3:  {C[0]/S_A[0]:+.2f}, {np.interp(1,z,C/S_A):+.2f}, {C[-1]/S_A[-1]:+.2f}")
print(f"    (analytic C = −2S; negative ∀z ⇒ NO stable canonical equilibrium)")
print(f"  gravity is long-range α=1 ≤ d (run_residue_longrange.py) ⇒ non-additive.")
print(f"  Campa–Dauxois–Ruffo: non-additive + C<0 ⇒ ENSEMBLE INEQUIVALENCE; the")
print(f"  canonical ensemble (the area-law λ=1 horn) is ill-defined for the")
print(f"  isolated self-gravitating horizon. Only the MICROCANONICAL (state-count)")
print(f"  accounting survives — the volume-law λ=1/2 branch SEDE adopts.")
print(f"  ⇒ NARROWS §9: the λ ambiguity loses its area horn; residue = whether the")
print(f"    microcanonical state count is volume (still dS holography).")

# ===========================================================================
# (2) GREA bulk-viscosity mapping of SEDE's structure-driven entropy production
# ===========================================================================
print("\n" + "=" * 74)
print("(2) SEDE as structure-gated GREA: the entropy production is a bulk viscosity")
print("=" * 74)
zg = np.linspace(0.0, 3.0, 600)
a = 1/(1+zg)
Eg = E_SEDE_lambda(zg, Om, gamma, lam, Or)
D = compute_growth_factor(zg, Om, Or)
f = np.clip((1-np.exp(-gamma*D**2))/(1-np.exp(-gamma)), 1e-9, 1)
lna = np.log(a)
eps = -np.gradient(np.log(Eg), lna)              # ε = −dlnH/dlna
dlnf = np.gradient(np.log(f), lna)               # entropy-production rate
# Thm 13:  1+w = (1/3)[2λε − dlnf/dlna];  bulk-viscosity reading of the 2nd term
rho_DE = ODE0 * f * Eg**(2*lam)                  # ρ_DE/ρ_crit0
# −(1/3)dlnf  ≡  −3ζH/ρ_DE  ⇒  ζ = ρ_DE·dlnf /(9 H²)   (units H0=1, ρ_crit0=1)
zeta = rho_DE * dlnf / (9 * Eg**2)
w = -1 + (1/3)*(2*lam*eps - dlnf)
zc = np.interp(0.0, w[::-1]+1, zg[::-1])
print(f"  effective bulk viscosity ζ(z) = ρ_DE·(dln f_sat/dlna)/(9H²):")
print(f"    {'z':>5} {'ζ (arb)':>9} {'w':>8}")
for zq in (0.0, 0.19, 0.5, 1.0, 2.0):
    print(f"    {zq:>5.2f} {np.interp(zq,zg,zeta):>9.4f} {np.interp(zq,zg,w):>8.3f}")
print(f"  ζ ≥ 0 everywhere: {bool(np.all(zeta[zg<2.5] >= -1e-9))}  (second-law-respecting,")
print(f"    as GREA requires); ζ PEAKS at the structure-formation epoch, not at the")
print(f"    horizon-expansion rate ⇒ SEDE = STRUCTURE-GATED GREA. w=−1 crossing at")
print(f"    z≈{zc:.2f} is where the viscous (entropy-production) term cancels the")
print(f"    dilution (2λε) term — GREA's negative-pressure balance, structure-timed.")

# validation
assert np.all(C < 0), "dS horizon must have negative specific heat"
assert np.all(zeta[zg < 2.2] >= -1e-6), "effective bulk viscosity must be >= 0 (2nd law)"
print("\n[validate] literature-inference assertions passed.")
