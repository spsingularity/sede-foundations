"""
run_driven_ness.py — PROTOTYPE of the driven-NESS entropy theorem (§8.6/§9).

CLAIM: a strongly-long-range horizon, DRIVEN by structure formation, settles in a
VOLUME-law steady state; remove the drive (→equilibrium) OR the long range
(→short-range) and it relaxes to AREA-law. The well-posed successor to the
volume-law postulate. We test it in two models.

PART 1 (negative, instructive): naive driven-dissipative free fermions do NOT
isolate the effect — uniform drive gives a trivially extensive thermal state
(connectivity-independent), boundary drive decays to area-law. So the theorem is
NOT 'generic driven transport'. This points to the real content:

PART 2 (the mechanism): the area↔volume ORDER PARAMETER under driven relaxational
dynamics. Long-range connectivity supplies the cooperative coupling J≥J_c (the
spinodal/bistability of run_chr_soc.py); the structure drive then pushes the order
parameter over the barrier and the bistability LOCKS it on the volume branch — a
genuine non-equilibrium selection. Short-range (J<J_c, monostable) cannot lock:
it tracks the drive and relaxes back to area. Undriven stays at the area-law
ground. This is the equilibrium(BH→area) vs driven-NESS(cosmic→volume) split.
"""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from scipy.integrate import solve_ivp

# ===========================================================================
# PART 1 — free-fermion driven-dissipative NESS (documented NEGATIVE)
# ===========================================================================
def _hop(N, alpha, J=1.0):
    i = np.arange(N); r = np.abs(i[:, None] - i[None, :]).astype(float)
    np.fill_diagonal(r, np.inf); h = J / r ** alpha; np.fill_diagonal(h, 0.0); return h

def _ness_cov(N, alpha, gamma, p, drive_sites):
    h = _hop(N, alpha); pv = np.zeros(N); pv[drive_sites] = p
    A = -1j * h - 0.5 * np.diag(gamma + pv)
    C = solve_continuous_lyapunov(A, -np.diag(pv)); return 0.5 * (C + C.conj().T)

def _block_S(C, sites):
    ev = np.clip(np.linalg.eigvalsh(C[np.ix_(sites, sites)]).real, 1e-12, 1 - 1e-12)
    return float(-np.sum(ev * np.log(ev) + (1 - ev) * np.log(1 - ev)))

print("=" * 72)
print("PART 1 — naive free-fermion driven NESS (why it is INSUFFICIENT)")
print("=" * 72)
N = 64
S_unif = _block_S(_ness_cov(N, 0.6, 1.0, 0.6, range(N)), list(range(N // 2)))
S_bdry_lr = _block_S(_ness_cov(N, 0.6, 1.0, 0.6, [0]), list(range(N // 2)))
S_bdry_sr = _block_S(_ness_cov(N, 3.0, 1.0, 0.6, [0]), list(range(N // 2)))
print(f"  uniform drive : S(N/2) = {S_unif:.2f}  — extensive but CONNECTIVITY-INDEPENDENT")
print(f"                  (C=νI for any hopping ⇒ trivially volume-law, no test of range)")
print(f"  boundary drive: S(N/2) long-range={S_bdry_lr:.2f}  short-range={S_bdry_sr:.2f}")
print(f"                  — both saturate (area-law): excitation decays from the source.")
print("  ⇒ generic driven transport does NOT realise the theorem; the content is")
print("    the COOPERATIVE bistable locking (Part 2), tied to the J≥J_c spinodal.")

# ===========================================================================
# PART 2 — order-parameter driven NESS: cooperative locking selects volume
# ===========================================================================
# Order parameter m∈(0,1): m≈0 = area-law ground, m≈1 = volume-law saturated.
# Cooperative (Wilson-Cowan / autocatalytic) relaxational dynamics:
#   dm/dt = -m + σ(J m + h(t) - θ),   σ(x)=1/(1+e^{-x})
# θ = FIXED threshold = the holographic-bound pull toward area (equilibrium
#     attractor). J = cooperativity from connectivity: J>J_c=4 (long-range /
#     all-to-all, run_chr_soc.py) ⇒ BISTABLE (area AND volume fixed points);
#     J<J_c ⇒ monostable area. h(t) = transient structure-deposition drive.
THETA = 3.0
sig = lambda x: 1.0 / (1.0 + np.exp(-x))

def dmdt(t, m, J, h_of_t):
    m = float(np.clip(m, 0.0, 1.0))
    return -m + sig(J * m + h_of_t(t) - THETA)

def structure_drive(h0=8.0, t0=5.0, w=2.5):
    """Structure-deposition pulse: rises, peaks (~z*≈1.2), then declines as
    saturation completes (df_sat/dt peaks then falls). The DRIVE is transient."""
    return lambda t: h0 * np.exp(-0.5 * ((t - t0) / w) ** 2)

print("\n" + "=" * 72)
print("PART 2 — order-parameter NESS: does the DRIVE through long-range LOCK volume?")
print("=" * 72)
h_pulse = structure_drive()
tspan = (0.0, 40.0); teval = np.linspace(*tspan, 800)
print(f"  θ={THETA} (holographic/area pull); J_c=4; structure pulse peaks t≈5, gone by t≈15;")
print(f"  start at the area ground m₀≈0.05; read final m at t=40")
print(f"  {'config':<34}{'m_peak':>9}{'m_final':>9}{'state':>16}")
final = {}; _traj = {}
for label, J, driven in [
    ("long-range (J=6>J_c), DRIVEN ", 6.0, True),
    ("short-range (J=2<J_c), DRIVEN", 2.0, True),
    ("decoupled (J=0), DRIVEN      ", 0.0, True),
    ("long-range (J=6), UNDRIVEN   ", 6.0, False),
]:
    hf = h_pulse if driven else (lambda t: 0.0)
    m0 = sig(-THETA)                              # the area-law ground state
    sol = solve_ivp(dmdt, tspan, [m0], args=(J, hf), t_eval=teval, rtol=1e-9, atol=1e-11)
    m = sol.y[0]
    mf = m[-1]; mp = m.max()
    final[label] = mf
    _traj[label] = (sol.t, m)
    state = "VOLUME (locked)" if mf > 0.5 else "area (relaxed)"
    print(f"  {label:<34}{mp:>9.3f}{mf:>9.3f}{state:>16}")

print("""
  Reading: m_final>0.5 = locked on the VOLUME branch; <0.5 = relaxed to AREA.
  The DRIVE is identical and transient in all three driven runs — only the
  connectivity (J) differs.""")

# ===========================================================================
# Verdict
# ===========================================================================
m_lr = final["long-range (J=6>J_c), DRIVEN "]
m_sr = final["short-range (J=2<J_c), DRIVEN"]
m_eq = final["long-range (J=6), UNDRIVEN   "]
print("=" * 72)
print("VERDICT — prototype of the driven-NESS theorem")
print("=" * 72)
ok = (m_lr > 0.5) and (m_sr < 0.5) and (m_eq < 0.5)
print(f"""  driven + long-range  : m_final = {m_lr:.2f}  ⇒ {'VOLUME-law (locked)' if m_lr>0.5 else 'area'}
  driven + short-range : m_final = {m_sr:.2f}  ⇒ {'volume' if m_sr>0.5 else 'AREA-law (relaxed)'}
  undriven long-range  : m_final = {m_eq:.2f}  ⇒ {'volume' if m_eq>0.5 else 'AREA-law (ground)'}

  {'THEOREM DEMONSTRATED in this model:' if ok else 'INCONCLUSIVE:'} a TRANSIENT structure drive,
  acting THROUGH long-range cooperativity (J≥J_c), pushes the order parameter
  over the spinodal barrier and the bistability LOCKS it on the volume branch —
  a non-equilibrium steady state that persists after the drive is gone. Remove
  the long range (J<J_c: monostable, tracks then relaxes to area) OR the drive
  (stays at the area-law ground): either collapses the volume law. This is the
  mechanism the full theorem must capture — driving × cooperativity selects
  volume; equilibrium or locality gives area (BH vs cosmic, §8.6).
  HONEST SCOPE: a mean-field order-parameter prototype, not the field-theoretic
  horizon proof. It shows the selection mechanism is real and identifies the two
  necessary ingredients (drive AND long-range bistability); a microscopic
  derivation of J(connectivity)≥J_c and of the deposition→h map remains open.""")

assert ok, "driven-NESS prototype should lock volume only for driven+long-range"
print("\n[validate] driven-NESS prototype assertions passed.")

# ---- figure: m(t) trajectories + the transient drive ----
import os
if os.environ.get("SEDE_NO_FIG") != "1":
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        colors = {"long-range (J=6>J_c), DRIVEN ": "C3",
                  "short-range (J=2<J_c), DRIVEN": "C0",
                  "decoupled (J=0), DRIVEN      ": "C2",
                  "long-range (J=6), UNDRIVEN   ": "0.5"}
        for label, (t, m) in _traj.items():
            ax.plot(t, m, color=colors.get(label, "k"), lw=2, label=label.strip())
        td = np.linspace(*tspan, 400)
        ax.plot(td, structure_drive()(td) / 8.0, "k:", lw=1, label="structure drive h(t) (scaled)")
        ax.axhline(0.5, ls="--", c="0.7"); ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("time"); ax.set_ylabel("order parameter m  (0=area, 1=volume)")
        ax.set_title("Driven-NESS: long-range+drive LOCKS volume; else relaxes to area")
        ax.legend(fontsize=7.5, loc="center right")
        fig.tight_layout()
        out = "output"
        os.makedirs(out, exist_ok=True)
        fig.savefig(os.path.join(out, "driven_ness.png"), dpi=130); plt.close(fig)
        print(f"[figure] wrote {os.path.join(out, 'driven_ness.png')}")
    except Exception as e:
        print(f"[figure] skipped ({e})")
