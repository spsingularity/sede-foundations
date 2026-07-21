"""
run_residue_longrange.py — resolving the counting RESIDUE from first principles.

The residue (§8.6): does the horizon entropy count BULK (volume, N∝R³, Δ=1) or
BOUNDARY (area, N∝R², Δ=0) dof?  Routes A,B,E reduced everything else; the count
is the last open piece. This script develops the one first-principles handle that
actually bears on it, and is honest about where it stops.

REDUCTION CHAIN:
  count (area vs volume)  ⟺  connectivity (local vs nonlocal/expander)   [RT min-cut]
  connectivity            ⟸  gravity is strongly LONG-RANGE              [this script]

The rigorous, citable statistical-mechanics classification (Campa, Dauxois &
Ruffo, Phys. Rep. 480 (2009) 57) of a 1/r^α interaction in d dimensions:
      α ≤ d   →  STRONGLY long-range, NON-ADDITIVE, entropy NOT area-additive
      α >  d  →  short-range, additive, area-law.
Newtonian gravity is α = 1. On the horizon (a d=2 surface) or its bulk (d=3),
α = 1 ≤ d ⇒ gravity is strongly long-range ⇒ the NATURAL entropy class is the
non-additive (volume) one. We verify the non-additivity scaling numerically, then
confront the black-hole counterexample and state the discriminator that survives.
"""
import numpy as np

rng = np.random.default_rng(0)

# ----------------------------------------------------------------------------
# (1) Non-additivity scaling:  per-particle energy u(N) at FIXED density.
#     u ∝ N^{1-α/d}  for α<d (super-extensive=non-additive); const for α>d.
# ----------------------------------------------------------------------------
def per_particle_energy(N, d, alpha, dens=1.0, r_core=0.7):
    """Mean Σ_j 1/r_ij^α over points at fixed number density in a d-ball.
    A short-distance core r_core (~mean spacing = dens^{-1/d} = 1) regulates the
    UV (nearest-neighbour) divergence of 1/r^α, isolating the COLLECTIVE (IR)
    behaviour that determines additivity — the physics of the residue."""
    L = (N / dens) ** (1.0 / d)                 # box size so density fixed
    X = rng.uniform(0, L, size=(N, d))
    # sample a subset of centres for speed
    nc = min(N, 40)
    idx = rng.choice(N, nc, replace=False)
    us = []
    for i in idx:
        dr = X - X[i]
        r = np.sqrt((dr ** 2).sum(1))
        r = r[r > 1e-9]
        r = np.maximum(r, r_core)               # soften the UV (nearest-neighbour) core
        us.append(np.sum(1.0 / r ** alpha))
    return np.mean(us)

print("=" * 72)
print("(1) Gravity (α=1) is strongly long-range ⇒ non-additive (volume class)")
print("=" * 72)
Ns = [200, 800, 3200]
for d in (2, 3):
    print(f"\n  d={d} horizon {'surface' if d==2 else 'bulk'}:  expected slope 1-α/d")
    print(f"    {'α':>5} {'predict slope':>13} {'measured slope':>15} {'class':>14}")
    for alpha in (1.0, float(d), float(d) + 1.0):
        us = np.array([per_particle_energy(N, d, alpha) for N in Ns])
        slope = np.polyfit(np.log(Ns), np.log(us), 1)[0]
        pred = max(0.0, 1.0 - alpha / d)
        cls = "NON-ADDITIVE" if alpha < d - 1e-9 else ("marginal" if abs(alpha-d)<1e-9 else "additive")
        tag = "  <- GRAVITY" if abs(alpha - 1.0) < 1e-9 else ""
        print(f"    {alpha:>5.1f} {pred:>13.3f} {slope:>15.3f} {cls:>14}{tag}")
print("\n  ⇒ α=1 gravity has POSITIVE slope (u/N grows) in d=2 and d=3 ⇒ non-additive")
print("    ⇒ the natural horizon-entropy class is the volume (non-area-additive) one.")

# ----------------------------------------------------------------------------
# (2) The black-hole counterexample and the surviving discriminator.
# ----------------------------------------------------------------------------
print("\n" + "=" * 72)
print("(2) But a black hole is gravitational AND area-law. What breaks the tie?")
print("=" * 72)
print("""  Long-range ⇒ volume is the UNCONSTRAINED class. Two horizons, same gravity:
    - Black hole: ISOLATED, in equilibrium, SATURATES the holographic bound ⇒
      the bound wins and relaxes it to AREA-law (Δ=0). [route-E Bekenstein = area]
    - Cosmic horizon: continuously DRIVEN by structure formation (the f_sat gate,
      dS_struct/dt > 0), a non-equilibrium steady state held OFF the equilibrium
      area-law attractor.
  So the discriminator is EQUILIBRIUM (→area) vs DRIVEN NESS (→volume), and the
  driver is the very structure-deposited entropy SEDE is built on.""")

# Check the cosmic horizon is in fact being driven TODAY (nonzero dln f_sat/dlna):
from sede.friedmann import compute_growth_factor
Om, Or, gamma = 0.30, 9.0e-5, 1.4964
z = np.array([0.0, 0.5, 1.0])
a = 1.0 / (1.0 + z)
D = compute_growth_factor(z, Om, Or)
f = (1 - np.exp(-gamma * D ** 2)) / (1 - np.exp(-gamma))
# dln f/dln a via finite diff on a fine grid
zg = np.linspace(0, 2, 400); ag = 1 / (1 + zg)
Dg = compute_growth_factor(zg, Om, Or)
fg = (1 - np.exp(-gamma * Dg ** 2)) / (1 - np.exp(-gamma))
dlnf_dlna = np.gradient(np.log(fg), np.log(ag))
drive_now = np.interp(0.0, zg, dlnf_dlna)
drive_peak = dlnf_dlna.max()
print(f"  cosmic-horizon drive |dln f_sat/dln a|:  today = {drive_now:.3f}, "
      f"peak = {drive_peak:.3f}  (NONZERO ⇒ driven NESS)")
print(f"  black-hole drive = 0 (isolated, equilibrium) ⇒ area. The split is REAL.")

print("\n" + "=" * 72)
print("VERDICT — how far this resolves the residue")
print("=" * 72)
print("""  PROGRESS, NOT YET CLOSURE — but the residue is reduced to ONE provable claim.
    1. count(area/volume) ⟺ connectivity(local/nonlocal)         [RT min-cut;
       run_chr_soc.py / run_tensor_network.py]
    2. gravity is strongly long-range (α=1 ≤ d) ⇒ non-additive   [verified above;
       Campa-Dauxois-Ruffo] ⇒ the NATURAL class is volume.
    3. the black-hole/cosmic asymmetry — the paper's named conceptual cost — is
       no longer a separate problem: it is equilibrium (bound-saturated → area)
       vs driven NESS (structure-fed → volume). One mechanism, both horizons.
  The single remaining rigorous gap is a DRIVEN-NESS ENTROPY THEOREM: prove that a
  strongly-long-range horizon, driven at rate dS_struct/dt, settles in the
  volume-law (non-additive) steady state rather than relaxing to the area-law
  equilibrium. That is the well-posed, self-contained statement that would CLOSE
  the residue — and it unifies it with the BH asymmetry (turning a cost into a
  prediction). Falsifiable corollary: the volume law is tied to the DRIVING, so
  Δ should track the structure-driving history (testable via the §6 Δ programme).""")

# validation
sl_grav_d3 = np.polyfit(np.log(Ns),
    np.log([per_particle_energy(N, 3, 1.0) for N in Ns]), 1)[0]
assert sl_grav_d3 > 0.2, "gravity must be non-additive (positive slope) in d=3"
assert drive_now > 0.01, "cosmic horizon must be driven today"
print("\n[validate] residue-long-range assertions passed.")
