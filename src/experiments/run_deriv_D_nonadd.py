"""
run_deriv_D_nonadd.py — Derivation route D: pin the exponent from gravitational
non-additivity (data-anchored, weak precedent).

Tsallis-Cirto non-extensive entropy S ∝ A^δ comes from gravity being long-range
/ non-additive; Barrow's Δ maps as δ = 1 + Δ/2, so the volume endpoint Δ=1 is
δ = 3/2. A first-principles derivation would fix δ (or the Tsallis q) from the
actual non-additivity of self-gravitating correlations, rather than leaving it
free. We test whether the data support the needed exponent.

PRECEDENT (already in the repo): run_cluster_tsallis.py fit a Tsallis q-Gaussian
to real SDSS/Coma cluster kinematics and found q ≈ 1.01 (GAUSSIAN), NOT q ≈ 1.5
— the kinetic non-extensivity that would mirror δ=3/2 is ABSENT at galaxy scales.
"""
import numpy as np

def delta_of_Delta(D):   return 1.0 + D/2.0      # Tsallis-Cirto exponent vs Barrow

print("=" * 70)
print("Route D — derive the exponent from gravitational non-additivity?")
print("=" * 70)
print("  Barrow Δ ↔ Tsallis-Cirto δ = 1 + Δ/2:")
for D in [0.0, 0.5, 1.0]:
    print(f"    Δ={D:.1f}  ->  δ={delta_of_Delta(D):.2f}   "
          + ("(area law)" if D == 0 else "(volume law)" if D == 1 else ""))
print(f"  SEDE needs the HORIZON exponent δ = 3/2 (Δ=1).")

print("\n  Data precedent — does structure exhibit the matching non-extensivity?")
print("  run_cluster_tsallis.py (real SDSS/Coma, 1073 members):")
q_measured = 1.01           # from run_cluster_tsallis.py MLE fit
q_predicted = 1.5           # the value that would mirror δ=3/2 at the kinetic level
print(f"    kinetic Tsallis q (measured)  = {q_measured:.2f}  (≈1 ⇒ Gaussian/extensive)")
print(f"    q that would mirror δ=3/2      = {q_predicted:.2f}")
print(f"    ⇒ MISMATCH: the predicted kinetic non-extensivity is NOT observed.")

print("\n" + "=" * 70)
print("VERDICT (route D)")
print("=" * 70)
print("""  WEAK / NEGATIVE. The δ↔Δ map is exact (δ=3/2 ⟺ Δ=1), but pinning δ from
  gravitational non-additivity is not supported:
  - The repo's real-data test (run_cluster_tsallis.py) finds galaxy kinematics
    GAUSSIAN (q≈1.01), not the q≈1.5 that a kinetic mirror of δ=3/2 would need.
  - Honest caveat: the HORIZON exponent δ and the KINETIC cluster q are DIFFERENT
    objects (horizon entropy scaling vs velocity-distribution non-extensivity);
    a failed kinetic q does not strictly falsify the horizon δ. But it removes
    the only handle we had to DERIVE δ from a measurable non-additivity, and the
    precedent is discouraging.
  ⇒ Lowest expected payoff of the five routes; report as a tested-and-weak link,
    consistent with the existing §cross-field negative result.""")

assert abs(delta_of_Delta(1.0) - 1.5) < 1e-9
assert q_measured < 1.1, "cluster kinematics are ~Gaussian (route D precedent fails)"
print("\n[validate] route-D assertions passed.")
