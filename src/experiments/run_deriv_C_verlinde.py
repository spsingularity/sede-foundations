"""
run_deriv_C_verlinde.py — Derivation route C: emergent-gravity unification.

Instead of postulating volume-counting for SEDE's dark ENERGY, postulate
entropic/emergent gravity (Verlinde 2016), which ALREADY invokes a volume-law
de Sitter entanglement entropy to produce apparent dark MATTER. If the SAME
volume-law entropy delivers both, the postulate is not new to SEDE — it is the
shared root of two independent dark-sector anomalies.

TEST: does one volume-law scale give BOTH ρ_DE ~ ρ_crit (SEDE) AND the MOND/RAR
acceleration a₀ ~ cH₀/2π (Verlinde DM)?  An O(1) match motivates (does not
derive) the postulate by tying it to an independent observable.
"""
import numpy as np

G    = 6.674e-11
c    = 2.99792458e8
H0   = 67.4 * 1000 / 3.0857e22          # s^-1
TWO_PI = 2 * np.pi

a_H0  = c * H0                          # Hubble acceleration
a0_obs = 1.2e-10                        # observed RAR/MOND scale (McGaugh-Lelli-Schombert 2016)

print("=" * 70)
print("Route C — one volume-law entropy for dark ENERGY (SEDE) and dark MATTER (Verlinde)?")
print("=" * 70)
print(f"  cH₀                 = {a_H0:.3e} m/s²")
print(f"  cH₀ / 2π            = {a_H0/TWO_PI:.3e} m/s²")
print(f"  observed a₀ (RAR)   = {a0_obs:.3e} m/s²")
print(f"  (cH₀/2π) / a₀_obs   = {a_H0/TWO_PI/a0_obs:.2f}   (O(1) ⇒ same scale)")

# SEDE side: the volume-law horizon entropy fixes ρ_DE ~ ρ_crit (Ω_DE ~ O(1), CKN).
# Verlinde side: the SAME de Sitter volume entropy S_DE ∝ V gives an entropic
# force with characteristic acceleration a₀ = cH₀/(2π) (his eq. for the dark
# elastic response). Both are the cosmological-horizon volume-law scale.
rho_crit = 3 * H0**2 / (8 * np.pi * G)            # kg/m^3
rho_DE   = 0.685 * rho_crit
print(f"\n  ρ_crit             = {rho_crit:.3e} kg/m³")
print(f"  ρ_DE (Ω_DE=0.685)  = {rho_DE:.3e} kg/m³   (SEDE: volume-law horizon entropy ⇒ Ω_DE~O(1))")
print(f"  both the DE density and a₀ trace the SAME cosmological-horizon volume scale.")

print("\n" + "=" * 70)
print("VERDICT (route C)")
print("=" * 70)
print(f"""  UNIFICATION MOTIVATION (not a derivation). One volume-law de Sitter
  entanglement entropy underlies BOTH SEDE's dark energy (ρ_DE~ρ_crit via CKN)
  AND Verlinde's apparent dark matter (a₀=cH₀/2π, here {a_H0/TWO_PI/a0_obs:.2f}×a₀_obs). If
  emergent gravity with a volume-law entropy is adopted as the deeper postulate,
  SEDE's dark-energy volume law is a COROLLARY, and the dark sector is unified.
  - Buys: ties the postulate to an independent anomaly (galaxy rotation/RAR);
    borrows an existing framework's motivation.
  - Risk: Verlinde's programme is itself contested and not a complete derivation;
    you inherit its open issues (covariant formulation, cluster-scale tension).
  ⇒ A motivation/unification route, strengthening §7 (Verlinde link) and §8.1;
    not a closure of the postulate.""")

assert 0.5 < a_H0/TWO_PI/a0_obs < 1.5, "cH0/2pi should be O(1) times a0"
print("\n[validate] route-C assertions passed.")
