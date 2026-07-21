"""
run_deriv_E_nogo.py — Derivation route E: the no-go lemma.

CLAIM TO TEST: can the volume-law postulate be derived from standard entropy
bounds (maximum-entropy + Bekenstein)?  ANSWER: NO — energy bounds give the
AREA law. This is a negative result that SHARPENS the postulate: it proves the
volume law needs a state/counting input beyond any energy bound, and pins the
size of the missing ingredient (a factor R/ℓ_P, the §8.3 'seam').

Physics (all analytic, no fit):
  cosmic apparent horizon, R_H = c/H,  E_hor = ρ_crit · V_H = c³/(2GH)
  Bekenstein:   S ≤ 2π R E /(ħc)  =  π R²/ℓ_P²  =  ¼ (A/ℓ_P²)   →  AREA LAW
  Volume law:   S_vol ~ V_H/ℓ_P³  =  (A/ℓ_P²)·(R/ℓ_P)/3·...     →  exceeds it by R/ℓ_P
"""
import numpy as np

G    = 6.674e-11          # m^3 kg^-1 s^-2
c    = 2.99792458e8       # m/s
hbar = 1.054571817e-34    # J s
H0   = 67.4 * 1000 / 3.0857e22   # s^-1  (67.4 km/s/Mpc)

lP = np.sqrt(hbar * G / c**3)     # Planck length
R  = c / H0                       # Hubble radius
A  = 4 * np.pi * R**2             # horizon area
V  = (4.0/3.0) * np.pi * R**3     # horizon volume

# Gibbons-Hawking (area) entropy in Planck units
S_area = A / (4 * lP**2)

# Bekenstein bound on the enclosed energy E_hor = rho_crit * V = c^3/(2GH)
rho_crit = 3 * H0**2 / (8 * np.pi * G) * c**2 / c**2   # mass density
E_hor = rho_crit * V * c**2                            # energy (J), = c^4/(2 G H0)? check below
E_hor_closed = c**5 / (2 * G * H0)                     # ρ_crit V = c³/(2GH) [mass] ⇒ E = c⁵/(2GH)
S_bek = 2 * np.pi * R * E_hor / (hbar * c)             # Bekenstein bound (nats)

# Volume-law (Planckian entropy density) — what SEDE postulates
S_vol = V / lP**3

print("=" * 70)
print("Route E — can max-entropy/Bekenstein DERIVE the volume law?  (no-go)")
print("=" * 70)
print(f"  Hubble radius R_H      = {R:.3e} m   = {R/lP:.3e} ℓ_P")
print(f"  E_hor (= ρ_crit V_H)   = {E_hor:.3e} J   (closed c⁴/2GH₀ = {E_hor_closed:.3e} J)")
print()
print(f"  S_area (Gibbons-Hawking) = A/4ℓ_P²   = {S_area:.3e}")
print(f"  S_Bekenstein = 2πRE/ħc               = {S_bek:.3e}")
print(f"  ratio S_Bek/S_area                   = {S_bek/S_area:.3f}   (≈1 ⇒ Bekenstein = AREA law)")
print()
print(f"  S_volume (Planckian density V/ℓ_P³)  = {S_vol:.3e}")
print(f"  ratio S_vol/S_area                   = {S_vol/S_area:.3e}")
print(f"  predicted excess R/ℓ_P  (the 'seam') = {R/lP:.3e}")
print(f"  match (S_vol/S_area)/(R/ℓ_P) = {(S_vol/S_area)/(R/lP):.3f}  (O(1) ⇒ excess IS R/ℓ_P)")

print("\n" + "=" * 70)
print("VERDICT (route E)")
print("=" * 70)
print(f"""  The Bekenstein bound on the horizon's ENERGY gives S ≤ ¼(A/ℓ_P²) = the
  AREA law (ratio {S_bek/S_area:.2f}). So no energy/maximum-entropy argument yields the
  volume law: it would VIOLATE the Bekenstein bound by a factor R/ℓ_P ~ {R/lP:.1e}.

  This is the §8.3 'seam' derived from first principles: the volume law is NOT
  reducible to an energy bound — it requires an independent STATE/COUNTING input
  (a Planckian-density, volume-law-entangled horizon). The no-go does not derive
  the postulate; it PROVES the postulate is irreducible to standard bounds and
  fixes the exact size (R/ℓ_P = √[ρ_crit/ρ_Planck]⁻¹) of the missing ingredient.
  ⇒ Add as a lemma to §8.3 (strengthens the honesty; costs nothing).""")

# validation
assert abs(S_bek / S_area - 1.0) < 0.6, "Bekenstein bound must be ~area-law"
# the seam: S_vol/S_area = (V/ℓ_P³)/(A/4ℓ_P²) = 4(V/A)/ℓ_P = (4/3)(R/ℓ_P), coefficient 4/3.
# (checkpoint review: the old assertion compared to 1/3 — the wrong coefficient — and was
# disabled with `or True`; this is now a real, passing test of the seam's R/ℓ_P scaling.)
assert abs((S_vol / S_area) / (R / lP) - 4.0/3.0) < 0.05, "S_vol/S_area must scale as (4/3)(R/ℓ_P)"
assert 1e60 < R/lP < 1e62, "R/ℓ_P should be ~1e61 (the seam)"
print("\n[validate] route-E assertions passed.")
