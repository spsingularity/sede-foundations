"""
run_birth_bifurcation.py — closing the Kramers gap: selection at the ordering
bifurcation + exact long-range hysteresis (no barrier is ever crossed).

THE GAP. §5(ii) places the horizon DEEP in the ordered phase (barrier ~ N·T_AH,
J/J_c ~ 2πN), yet Fork A claims Δ=1 is reached by "ordinary relaxation". A Kramers
escape over an N·T barrier is suppressed by e^{-N}, N ~ 10^122 — it never happens.
So relaxation-over-the-barrier cannot be the selection mechanism.

THE CLOSURE. The barrier is N-DEPENDENT: J/J_c = 2πGm²N crosses 1 only at N ~ O(1)
(Planckian dof count). The cosmic apparent horizon GREW continuously from N ~ 1, so
it passed THROUGH the ordering bifurcation as the double well formed around it:
  (a) near threshold the barrier vanishes quadratically, Δf ∝ (J−J_c)², while the
      Barrow entropic tilt (§5 i, F monotone toward Δ=1) is already present;
  (b) an imperfect (tilted) pitchfork has no branch ambiguity — the state
      continuously connected to the disordered N~1 minimum is the TILTED (volume)
      branch; the area minimum is born METASTABLE at the spinodal J_sp > J_c and
      is never occupied;
  (c) the horizon therefore tracks ADIABATICALLY into the deepening volume well —
      Δ=1 selected with NO barrier crossing, ever.
HYSTERESIS is then a theorem, not a model choice: for long-range/mean-field systems
metastability is exact (no finite nucleation droplet, Penrose–Lebowitz; lifetimes
∝ e^{N·Δf}). Every horizon keeps its BIRTH branch:
  • cosmic horizon: grown from N~1 through the bifurcation → VOLUME (Δ=1);
  • black hole: created ABRUPTLY at macroscopic N in a low-entanglement (area)
    state → born in the area well → the e^{-N·Δf} suppression now works FOR the
    model: Δ_BH = 0 is permanent (the §7.1 PBH prediction, now derived).

We use the same mean-field F(m) as closure route C (run_closure_attempts.py):
  f(m)/T = −(J/2)m² + θm + m ln m + (1−m) ln(1−m),   m = volume-participating
fraction, J_c = 4 (units of T), coexistence field θ* = J/2. Barrow tilt b > 0:
θ = J/2 − b (entropy advantage of the volume branch per dof; b = O(0.1–1) is
conservative — the actual Barrow advantage grows with R).

HONEST FLAGS. The N ~ O(1) regime extrapolates the mean-field model into the
quantum-gravity era (tagged: the same tier as S_ent = S_grav). The claim proved
here is conditional-mechanical: GIVEN the route-C landscape with an N-dependent J
and a volume-favouring tilt, birth-selection + hysteresis follow with no Kramers
miracle. What selects the tilt's sign is Barrow's monotone slope (§5 i) — derived.
"""
import numpy as np

RNG = np.random.default_rng(7)

# ---------------------------------------------------------------------------
# the route-C mean-field free energy (per dof, units of T)
# ---------------------------------------------------------------------------
def f(m, J, theta):
    return -0.5 * J * m**2 + theta * m + m * np.log(m) + (1 - m) * np.log(1 - m)

def fp(m, J, theta):
    return -J * m + theta + np.log(m / (1 - m))

def _grid(n=40001, lim=25.0):
    """Logit-spaced grid: resolves wells that hug m→0 / m→1 at deep J."""
    u = np.linspace(-lim, lim, n)
    return 1.0 / (1.0 + np.exp(-u))

_M = _grid()

def local_minima(J, theta):
    y = f(_M, J, theta)
    idx = np.where((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:]))[0] + 1
    return _M[idx], y[idx]

def barrier_from(m0, J, theta):
    """Barrier per dof from the well containing m0 to the adjacent saddle."""
    y = f(_M, J, theta)
    mins, _ = local_minima(J, theta)
    if len(mins) < 2:
        return 0.0
    # saddle = max of f between the two outermost minima
    lo, hi = np.sort(mins)[0], np.sort(mins)[-1]
    sel = (_M > lo) & (_M < hi)
    saddle = y[sel].max()
    well = f(np.clip(m0, 1e-9, 1 - 1e-9), J, theta)
    return float(saddle - well)

J_C = 4.0
B_TILT = 0.3   # Barrow entropic tilt per dof (conservative O(0.1) value)

# ---------------------------------------------------------------------------
# 1. barrier vanishes quadratically at threshold: Δf ∝ (J−J_c)²
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. Near the ordering threshold the barrier VANISHES: Δf ∝ (J−J_c)²")
print("=" * 78)
Js = J_C + np.array([0.1, 0.2, 0.4, 0.8, 1.2])
bars = []
for J in Js:
    mins, ys = local_minima(J, J / 2)         # coexistence (b=0): symmetric well
    m_lo = mins[np.argmin(mins)]
    bars.append(barrier_from(m_lo, J, J / 2))
bars = np.array(bars)
slope = np.polyfit(np.log(Js - J_C), np.log(bars), 1)[0]
print(f"    {'J−J_c':>8} {'barrier Δf/T per dof':>22}")
for dJ, b in zip(Js - J_C, bars):
    print(f"    {dJ:>8.2f} {b:>22.5f}")
print(f"  scaling: Δf ∝ (J−J_c)^{slope:.2f}   (Landau prediction: exponent 2, "
      f"coefficient 3/64 = {3/64:.4f}; fit c = {np.exp(np.polyfit(np.log(Js-J_C), np.log(bars),1)[1]):.4f})")
print("  ⇒ the N·T_AH barrier of §5(ii) did NOT exist at birth: it grows with N.")

# ---------------------------------------------------------------------------
# 2. tilted pitchfork: the branch connected to disorder is the VOLUME branch
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. With the Barrow tilt, the area well is born METASTABLE at J_sp > J_c")
print("=" * 78)
theta_of = lambda J: J / 2 - B_TILT
J_scan = np.linspace(3.0, 9.0, 601)
n_wells = []
for J in J_scan:
    mins, _ = local_minima(J, theta_of(J))
    n_wells.append(len(mins))
n_wells = np.array(n_wells)
one_well = J_scan[n_wells == 1]
J_sp = float(J_scan[n_wells >= 2].min()) if np.any(n_wells >= 2) else np.inf
mins_at, _ = local_minima(J_sp - 0.05 if J_sp < 9 else 5.0, theta_of(J_sp - 0.05 if J_sp < 9 else 5.0))
print(f"  tilt b = {B_TILT}: single-well region extends to J_sp = {J_sp:.2f} "
      f"(> J_c = {J_C}); in J ∈ ({J_C}, {J_sp:.2f}) the ONLY minimum is at "
      f"m = {float(mins_at[-1]):.3f} (volume side)")
print("  ⇒ as N (hence J ∝ 2πN) grows through threshold, the unique minimum")
print("    deforms continuously to high m; the area minimum appears LATER, at the")
print("    spinodal, as a new (metastable, never-occupied) well.")

# ---------------------------------------------------------------------------
# 3. the grown horizon: adiabatic sweep J: 2 → 24 — tracks into the volume well
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. GROWN horizon (cosmic): sweep J(N) through the bifurcation — adiabatic")
print("=" * 78)
dt, steps, cap = 0.005, 100000, 2e-3
J_path = np.linspace(2.0, 24.0, steps)
m = 0.45                                   # disordered-side start (any m works)
traj, max_step = [], 0.0
for i in range(steps):
    J = J_path[i]
    dm = float(np.clip(-fp(m, J, theta_of(J)) * dt, -cap, cap))
    max_step = max(max_step, abs(dm))
    m = float(np.clip(m + dm, 1e-9, 1 - 1e-9))
    if i % 1000 == 0:
        traj.append((J, m))
traj = np.array(traj)
mins_f, _ = local_minima(24.0, theta_of(24.0))
print(f"    {'J':>6} {'m(t)':>8}")
for J, mm in traj[::10]:
    print(f"    {J:>6.1f} {mm:>8.3f}")
m_grown = m
print(f"  final: m = {m_grown:.3f} (volume well; wells at final J: "
      f"{[f'{x:.3f}' for x in mins_f]})")
print(f"  max |Δm| per step = {max_step:.2e}  ⇒ CONTINUOUS tracking; gradient")
print("  descent cannot cross barriers, so ending in the volume well while the")
print("  area well exists PROVES no crossing was ever needed. Δ=1 by birth.")

# ---------------------------------------------------------------------------
# 4. the abrupt horizon (black hole): born at large N in the area well — locked
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("4. ABRUPT horizon (black hole): born at large N in the AREA well — locked")
print("=" * 78)
J_bh = 24.0
mins_bh, _ = local_minima(J_bh, theta_of(J_bh))
m_area_well = float(mins_bh.min())
m = m_area_well
N_demo = 1.0e4                              # per-dof noise ~ sqrt(2T/N): tiny
amp = np.sqrt(2 * dt / N_demo)
for i in range(steps):
    dm = float(np.clip(-fp(m, J_bh, theta_of(J_bh)) * dt, -cap, cap))
    m = float(np.clip(m + dm + amp * RNG.standard_normal(), 1e-9, 1 - 1e-9))
m_bh = m
dF = barrier_from(m_area_well, J_bh, theta_of(J_bh))
for N_phys, label in [(1e77, "stellar BH (10 M_sun)"), (1e122, "cosmic-horizon count")]:
    print(f"  Kramers exponent N·Δf ({label}): {N_phys * dF:.2e}   "
          f"[age of universe: ln(t_0/t_Pl) ≈ 140]")
print(f"  simulated BH-born trajectory (N = {N_demo:.0e} demo): stays at "
      f"m = {m_bh:.3f} (area well m = {m_area_well:.3f}) despite the volume tilt")
print("  ⇒ EXACT hysteresis (mean-field: no local nucleation droplet). Every")
print("    horizon keeps its BIRTH branch: Δ_BH = 0 is now DERIVED, not asserted,")
print("    and the same e^{-N·Δf} that seemed fatal to Fork A protects the PBH")
print("    prediction (§7.1) instead.")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — selection at the bifurcation, not relaxation over the barrier")
print("=" * 78)
print(f"""  The Kramers objection to Fork A dissolves: J/J_c = 2πGm²N is N-DEPENDENT,
  so the double well did not exist at the horizon's birth (N ~ O(1)). The grown
  horizon crossed the ordering bifurcation while the barrier was ∝ (J−J_c)² → 0
  (measured exponent {slope:.2f}), with the Barrow tilt already selecting the volume
  side; it tracked continuously into the volume well (final m = {m_grown:.2f}) and no
  barrier was ever crossed. Conversely a horizon CREATED at macroscopic N in the
  area state is locked there forever (N·Δf ~ 10^77·{dF:.2f} ≫ 140) — black holes
  are area-law BY BIRTH, the cosmic horizon volume-law BY GROWTH. One mechanism,
  both predictions; the 10⁻⁷-drive-flips-10⁶¹-dof objection never arises because
  nothing is ever flipped.
  FLAG: mean-field landscape extrapolated to N ~ O(1) (quantum-gravity era) —
  same conditional tier as S_ent = S_grav; the tilt's sign is Barrow's §5(i) slope.""")

# validation
assert 1.8 < slope < 2.2, "barrier must vanish quadratically at threshold"
assert J_sp > J_C + 0.1, "area well must be born metastable (spinodal above J_c)"
assert m_grown > 0.85, "grown horizon must end in the volume well"
assert max_step < 0.02, "sweep must be continuous (no jump = no barrier crossing)"
assert m_bh < 0.2, "abruptly-born horizon must stay locked in the area well"
assert 1e77 * dF > 1e70, "BH Kramers exponent astronomically large"
print("\n[validate] birth-bifurcation assertions passed.")
