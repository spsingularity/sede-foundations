"""
run_neff_gate.py — deriving the gate: f_sat's functional form from activation
kinetics, and why the equilibrium mean-field alternative is REJECTED.

THE GAP. §7.1's state-dependent count writes N_eff = N_area + (N_vol−N_area)·f_sat
and admits this is "an effective-model statement one layer out" — an extra
postulate introduced to reconcile the equilibrium-area results with the driven
volume horizon. Separately, Paper I's f_sat(z) = (1−e^{−γD²})/(1−e^{−γ}) is a
stated CONVENTION (sede/gamma_computation.py), not a derived form.

THE CLOSURE (two results):
  (1) KINETIC DERIVATION OF THE GATE. If each structure deposit (violent-
      relaxation event, §5) activates a fraction of the REMAINING inactive volume
      capacity proportional to the deposit — the minimal absorption kinetics
          df = γ (1 − f) d(D²),      f(D→0) = 0
      — then f(D) = 1 − e^{−γD²}, i.e. EXACTLY Paper I's gate up to its z=0
      normalisation. The logistic form is the survival function of capacity
      activation under cumulative deposition; γ is the activation cross-section
      per unit growth. One line of kinetics replaces a convention.
  (2) THE GATE IS KINETIC, NOT EQUILIBRIUM. The alternative — read f_sat as the
      equilibrium global-minimum activation m*(θ(z)) of the route-C landscape —
      FAILS: deep in the ordered phase (J ≫ J_c) the global minimum switches
      discontinuously at the coexistence drive, giving a step-like activation,
      not the smooth growth-tracking gate. The smoothness of f_sat is itself
      evidence that the horizon activation is a DRIVEN, cumulative (NESS)
      process — consistent with the paper's framing — and under the birth-
      selection picture (run_birth_bifurcation.py) the BRANCH is fixed at birth
      while f_sat is the OCCUPANCY of the volume capacity by deposited structure
      entropy. Equilibrium horizons (BH, eternal dS) have D ≡ 0 ⇒ f = 0 ⇒ the
      area baseline of §7.1, with no extra postulate.
ADIABATICITY: activation tracks D²(z) because deposition is prompt —
t_dyn/t_H = Δ_vir^{−1/2} ≈ 0.075 ≪ 1 (run_violent_relaxation.py).
"""
import numpy as np
from scipy.integrate import solve_ivp
from sede.friedmann import compute_growth_factor

GAMMA = 1.5
OMEGA_M = 0.3

# ---------------------------------------------------------------------------
# 1. the kinetic derivation:  df = γ(1−f) d(D²)  ⇒  f = 1 − e^{−γD²}
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. Gate from activation kinetics:  df = γ(1−f) d(D²)  ⇒  f = 1 − e^{−γD²}")
print("=" * 78)
z_grid = np.linspace(0.0, 30.0, 400)
D = compute_growth_factor(z_grid, OMEGA_M)
D2_of_z = lambda z: np.interp(z, z_grid, D**2)

def rhs(z, f):
    # d(D²)/dz via finite difference of the tabulated growth factor
    dz = 1e-3
    dD2 = (D2_of_z(z + dz) - D2_of_z(max(z - dz, 0.0))) / (dz + min(z, dz))
    return GAMMA * (1.0 - f) * dD2

sol = solve_ivp(rhs, [30.0, 0.0], [1.0 - np.exp(-GAMMA * D2_of_z(30.0))],
                t_eval=z_grid[::-1], rtol=1e-9, atol=1e-12)
f_ode = sol.y[0][::-1]                       # back to ascending z
f_closed = 1.0 - np.exp(-GAMMA * D**2)
err = np.max(np.abs(f_ode - f_closed))
# Paper-I normalised gate
f_sat = f_closed / (1.0 - np.exp(-GAMMA))
print(f"  ODE solution vs closed form 1 − e^(−γD²):  max |err| = {err:.2e}")
print(f"    {'z':>6} {'D(z)':>8} {'f (kinetic)':>12} {'f_sat (Paper I)':>16}")
for zz in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]:
    k = np.argmin(np.abs(z_grid - zz))
    print(f"    {zz:>6.1f} {D[k]:>8.3f} {f_ode[k]:>12.4f} {f_sat[k]:>16.4f}")
print("  ⇒ Paper I's gate (1−e^{−γD²})/(1−e^{−γ}) IS the solution of the minimal")
print("    absorption kinetics: each deposit activates γ·d(D²) of the remaining")
print("    inactive capacity. The 'convention' is now a derivation; γ = activation")
print("    cross-section per unit growth² (the Paper-I Theorem-4C weight).")

# ---------------------------------------------------------------------------
# 2. the equilibrium alternative is step-like ⇒ rejected
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. Equilibrium m*(θ(z)) is STEP-like — the gate cannot be an equilibrium map")
print("=" * 78)
def f_land(m, J, theta):
    return -0.5 * J * m**2 + theta * m + m * np.log(m) + (1 - m) * np.log(1 - m)

def m_global(J, theta, n=8001):
    m = np.linspace(1e-6, 1 - 1e-6, n)
    y = f_land(m, J, theta)
    return float(m[np.argmin(y)])

J_DEEP = 8.0
B_SCALE = 1.2                                 # drive scale: theta = J/2 − b·D²
D2_axis = np.linspace(0.0, 1.0, 300)
m_eq = np.array([m_global(J_DEEP, J_DEEP / 2 - B_SCALE * d2) for d2 in D2_axis])
act_eq = (m_eq - m_eq[0]) / (m_eq[-1] - m_eq[0])
act_kin = (1.0 - np.exp(-GAMMA * D2_axis)) / (1.0 - np.exp(-GAMMA))
slope_eq = np.max(np.abs(np.diff(act_eq))) / np.mean(np.diff(D2_axis))
slope_kin = np.max(np.abs(np.diff(act_kin))) / np.mean(np.diff(D2_axis))
print(f"  deep-ordered landscape (J = {J_DEEP} ≫ J_c = 4), drive θ = J/2 − b·D²:")
print(f"    max d(activation)/d(D²):  equilibrium = {slope_eq:.1f}   kinetic gate = {slope_kin:.2f}")
print(f"    equilibrium activation jumps {act_eq.max()-act_eq.min():.2f} across "
      f"ΔD² ≈ {np.mean(np.diff(D2_axis)):.3f} — a STEP at coexistence (first-order switch)")
print("  ⇒ a smooth, growth-tracking f_sat CANNOT be the equilibrium global-minimum")
print("    map of a deep-ordered landscape; it must be (and is) the CUMULATIVE")
print("    kinetic occupancy of the volume capacity — the NESS reading, derived.")

# ---------------------------------------------------------------------------
# 3. limits and adiabaticity
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. Limits: equilibrium horizons f→0 (area baseline); tracking is prompt")
print("=" * 78)
Delta_vir = 178.0
t_ratio = Delta_vir ** -0.5
print(f"  D → 0 (black hole, eternal dS: no growth): f = {1-np.exp(-GAMMA*0.0):.1f} ⇒ N_eff = N_area")
print(f"    — §7.1's equilibrium-area baseline with NO extra postulate.")
print(f"  D → 1 (today): f_sat = {f_sat[0]:.3f} ⇒ N_eff → N_vol (driven, locked).")
print(f"  adiabaticity: t_dyn/t_H = Δ_vir^(−1/2) = {t_ratio:.3f} ≪ 1 — deposits are")
print(f"    prompt on the Hubble time, so f tracks D²(z) without lag.")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — N_eff(f_sat) demoted from postulate to kinetics")
print("=" * 78)
print(f"""  (1) Paper I's gate form is DERIVED: df = γ(1−f)d(D²) ⇒ f = 1−e^(−γD²)
      (ODE vs closed form: {err:.1e}). The §7.1 interpolation N_eff = N_area +
      (N_vol−N_area)·f_sat is then the occupancy statement of the same kinetics,
      not an independent assumption; equilibrium horizons are its D=0 limit.
  (2) The equilibrium (global-minimum) reading of the gate is REJECTED by shape:
      step-like ({slope_eq:.0f} vs {slope_kin:.1f} max slope) — the horizon activation is a
      driven cumulative process, which is the NESS picture stated mechanically.
  FLAGS: 'capacity activated ∝ deposit' is the minimal (linear) kinetics — the
  honest content of the driven-NESS conjecture, now localised in one line; the
  branch itself is fixed at birth (run_birth_bifurcation.py), and f_sat is the
  occupancy of that branch, not the branch selector.""")

# validation
assert err < 1e-4, "kinetic ODE must reproduce the closed-form gate"
assert slope_eq > 5 * slope_kin, "equilibrium activation must be step-like vs kinetic"
assert abs(f_sat[0] - 1.0) < 1e-6, "gate normalised to 1 at z=0"
assert f_ode[-1] < 0.05, "gate must vanish at high z (D -> 0, area baseline)"
assert t_ratio < 0.1, "deposition must be prompt (adiabatic tracking)"
print("\n[validate] neff-gate assertions passed.")
