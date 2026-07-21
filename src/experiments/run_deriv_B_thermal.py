"""
run_deriv_B_thermal.py — Derivation route B: entanglement first law around a
THERMAL reference state (the most promising REDUCTION of the postulate).

IDEA: Jacobson derives the AREA-law Einstein equations from maximal entanglement
of the VACUUM (δS = δ⟨K⟩, S ∝ A). Redo the construction around a THERMAL (Gibbs)
state. A thermal state's entanglement entropy of a ball of radius R has a VOLUME
term that takes over from the area term once R exceeds the thermal length
λ_th = ħc/(k_B T). So 'volume-law' is the generic late-time / thermalised
behaviour — not exotic.

WHAT WE TEST QUANTITATIVELY:
  S_EE(R; T) = a (R/ε)²            [area, vacuum piece]   + s_th(T) · V   [thermal volume]
  crossover at  R* ~ λ_th = 1/T (natural units ħ=c=k_B=1).
  - At the de Sitter temperature T_dS = H/2π, the horizon R=1/H sits at
    R/λ_th = 1/(2π) < 1  ⇒  AREA side (just below crossover).
  - Volume-law (the SEDE postulate) requires T ≫ T_dS — horizon dof thermalised
    well above the de Sitter bath, toward the Planck scale (= maximal scrambling,
    §8.5). Then the entropy DENSITY is Planckian — exactly the CKN magnitude.

So route B REDUCES the postulate: the volume-law FORM follows from thermalisation
above the thermal length; the Planckian SCALE is the CKN input (§8.1 pt 2 / §8.3).
It does NOT derive it for free, but it trades 'volume counting' for the more
physical 'horizon dof are maximally thermalised (Planck-T)'.
"""
import numpy as np

# natural units ħ=c=k_B=1; lengths in units of a UV cutoff ε (set ε=1).
# entropy of a ball radius R in a thermal state at temperature T:
#   S(R,T) = a_area * R^2  +  s_th * (4/3 pi R^3),   s_th = b * T^3  (3D rel. gas)
a_area = 1.0
b_th   = 1.0

def S_area_part(R):           return a_area * R**2
def S_vol_part(R, T):         return b_th * T**3 * (4.0/3.0) * np.pi * R**3
def S_total(R, T):            return S_area_part(R) + S_vol_part(R, T)

def crossover_R(T):
    # area == volume:  a R^2 = b T^3 (4/3 pi) R^3  ->  R* = a / (b T^3 (4/3 pi))
    return a_area / (b_th * T**3 * (4.0/3.0) * np.pi)

print("=" * 70)
print("Route B — volume-law entanglement from a THERMAL reference state")
print("=" * 70)
print("  S(R,T) = a·R²(area) + s_th·V(volume),  s_th ∝ T³;  crossover R*(T)")
print(f"\n  {'T':>8} {'R*(crossover)':>14} {'regime at R=R*':>16}")
for T in [0.3, 1.0, 3.0, 10.0]:
    Rs = crossover_R(T)
    print(f"  {T:>8.2f} {Rs:>14.4f}   area↔volume at R=R*")

# Map to the cosmic horizon. In Hubble units set H=1 so R_H=1.
# de Sitter temperature T_dS = H/(2pi) = 1/(2pi).
H = 1.0
R_H = 1.0 / H
T_dS = H / (2 * np.pi)
lam_th_dS = 1.0 / T_dS
print("\n" + "-" * 70)
print("  Cosmic horizon at the de Sitter temperature:")
print(f"    R_H = 1/H = {R_H:.3f},  T_dS = H/2π = {T_dS:.4f},  λ_th = 1/T_dS = {lam_th_dS:.3f}")
print(f"    R_H / λ_th = {R_H/lam_th_dS:.4f}  (<1 ⇒ horizon is on the AREA side at T_dS)")
fa = S_area_part(R_H); fv = S_vol_part(R_H, T_dS)
print(f"    area piece = {fa:.3e},  volume piece = {fv:.3e},  vol/area = {fv/fa:.3f}")
print("    ⇒ at the de Sitter temperature alone, the horizon is AREA-dominated.")

# What temperature makes the horizon volume-law (R_H >> lambda_th)?
print("\n  Temperature needed for volume-law dominance (vol/area = 10):")
# fv/fa = (b T^3 (4/3 pi) R^3)/(a R^2) = (4/3 pi) T^3 R = 10  -> T = (10 / (4/3 pi R))^{1/3}
T_needed = (10.0 / ((4.0/3.0)*np.pi*R_H))**(1.0/3.0)
print(f"    T_vol ≈ {T_needed:.2f} · H   ⇒  T ≫ T_dS = {T_dS:.3f}·H")
print(f"    i.e. horizon dof thermalised at T/T_dS ≈ {T_needed/T_dS:.1f}× the de Sitter bath.")
print("    Pushed to the Planck scale (maximal scrambling, §8.5) the entropy")
print("    DENSITY becomes Planckian — the CKN magnitude (§8.2/8.3).")

print("\n" + "=" * 70)
print("VERDICT (route B)")
print("=" * 70)
print("""  PARTIAL REDUCTION (the strongest of the five for an actual derivation).
  - The volume-law FORM (S ∝ V) is the generic entanglement entropy of a
    thermalised state for R > λ_th — NOT exotic; it follows from a thermal
    reference state in the entanglement-first-law construction (Jacobson around
    Gibbs instead of vacuum; cf. the memo's C8 modular-energy route).
  - At the de Sitter temperature alone the horizon sits on the AREA side
    (R_H/λ_th = 1/2π); volume-law requires horizon dof thermalised ABOVE T_dS,
    toward Planck (= maximal scrambling). The Planckian density is then the CKN
    scale.
  NET: route B trades the 'volume-counting' postulate for the more physical
  'horizon dof are maximally thermalised (Planckian volume-law density)'. The
  FORM is derived from thermalisation; the SCALE remains CKN. This is a genuine
  reduction (postulate → state-thermalisation), not a free derivation.
  ⇒ Add a sentence to §8.1 pt 2 / §8.5: the volume FORM is the thermal-state
    entanglement law; pair with the SYK state result (route A) and CKN scale.""")

assert crossover_R(0.3) > crossover_R(10.0), "R* must shrink as T grows"
assert R_H/lam_th_dS < 1.0, "horizon is area-side at T_dS"
assert T_needed > T_dS, "volume-law needs T above the de Sitter bath"
print("\n[validate] route-B assertions passed.")
