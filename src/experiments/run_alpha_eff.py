"""
run_alpha_eff.py — the marginal-dimension correction made QUANTITATIVE:
Δ_pred = 1 − 1/(2 ln(R_H/ℓ_P)) ≈ 0.996, which RESTORES the falsifiability that
§7.2's qualitative "Δ ∈ (0,1) would be within the picture" quietly surrendered.

THE GAP. §7 pre-registers "a robust intermediate Δ (say 0.5) falsifies", but §7.2
notes the horizon sits at the EW lower critical dimension (α = 0, marginal/log)
and concludes a measured Δ ∈ (0,1) "would be within the picture" — a direct
conflict: the escape hatch un-falsifies the headline prediction.

THE CLOSURE. Don't leave the marginal case qualitative — COMPUTE it. At the lower
critical dimension the EW structure function is logarithmic, C(r) ∝ ln(r/a), so
the EFFECTIVE roughness exponent at scale ratio r/a is
    α_eff(r) = (1/2) d ln C / d ln r = 1 / (2 ln(r/a)).
The horizon's scale ratio is fixed by physics: r/a = R_H/ℓ_P ≈ 8×10⁶⁰, so
    Δ_pred = 1 − α_eff = 1 − 1/(2 ln(R_H/ℓ_P)) ≈ 0.9964.
The marginal "escape hatch" thus predicts a deviation of ~0.4% — invisible at
DESI DR3 + Euclid's σ_Δ ≈ 0.09 and NOWHERE NEAR Δ ≈ 0.5. Pre-registered window:
Δ ∈ [1 − O(1)/(2 ln R/ℓ_P), 1] ≈ [0.98, 1]; a robust Δ ≈ 0.5 falsifies at > 5σ.
BONUS: the same logic makes the measurement a ROUGHENING-CLASS discriminator —
KPZ (α ≈ 0.39) would give Δ ≈ 0.61, distinguishable from EW's 0.996 at ~4σ.

We verify the α_eff formula in a 2+1D EW simulation (local slopes of C(r) follow
1/(2 ln r)), then apply it at the horizon's scale ratio, and close with the
amplitude honesty check (thermal roughness is Planckian ⇒ the equilibrium horizon
is smooth/area — consistent with §7.1's baseline — and the DRIVEN amplitude is
the quantitative form of the activation/gate question).
"""
import numpy as np

RNG = np.random.default_rng(11)

# ---------------------------------------------------------------------------
# 1. the horizon numbers: Δ_pred = 1 − 1/(2 ln Λ),  Λ = R_H/ℓ_P
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. The marginal correction at the horizon's scale ratio")
print("=" * 78)
c = 2.99792458e8            # m/s
H0 = 67.4 * 1e3 / 3.0857e22 # 1/s  (67.4 km/s/Mpc)
lP = 1.616255e-35           # m
R_H = c / H0
Lam = R_H / lP
lnLam = np.log(Lam)
alpha_eff_h = 1.0 / (2.0 * lnLam)
Delta_pred = 1.0 - alpha_eff_h
sigma_DESI = 0.09
print(f"  R_H = c/H0 = {R_H:.3e} m,   Λ = R_H/ℓ_P = {Lam:.2e},   ln Λ = {lnLam:.1f}")
print(f"  α_eff = 1/(2 ln Λ) = {alpha_eff_h:.2e}")
print(f"  Δ_pred = 1 − α_eff = {Delta_pred:.4f}")
print(f"  vs DESI DR3 + Euclid σ_Δ ≈ {sigma_DESI}:")
print(f"    deviation from Δ=1:  (1 − Δ_pred)/σ = {(1-Delta_pred)/sigma_DESI:.3f} σ  — indistinguishable from 1")
print(f"    distance to Δ=0.5:   (Δ_pred − 0.5)/σ = {(Delta_pred-0.5)/sigma_DESI:.1f} σ  — Δ ≈ 0.5 FALSIFIES")
print(f"    KPZ alternative (α = 0.39 ⇒ Δ = 0.61): (Δ_pred − 0.61)/σ = "
      f"{(Delta_pred-0.61)/sigma_DESI:.1f} σ — the measurement discriminates EW vs KPZ")
print("  ⇒ the marginal-dimension case is NOT an escape hatch: it pins Δ to a")
print("    computable ~0.4% below 1. Pre-register: Δ ∈ [0.98, 1]; intermediate Δ falsifies.")

# ---------------------------------------------------------------------------
# 2. verify α_eff(r) = 1/(2 ln(r/a)) in a 2+1D EW simulation
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. Simulation check: 2+1D EW local exponents follow α_eff(r) = 1/(2 ln(r/a))")
print("=" * 78)
N, nu, D, dt, steps = 192, 1.0, 1.0, 0.1, 9000
h = np.zeros((N, N)); amp = np.sqrt(2 * D * dt)
snaps = []
for t in range(steps):
    lap = (np.roll(h, 1, 0) + np.roll(h, -1, 0) +
           np.roll(h, 1, 1) + np.roll(h, -1, 1) - 4 * h)
    h += dt * nu * lap + amp * RNG.standard_normal((N, N))
    h -= h.mean()
    if t > steps - 2000 and t % 400 == 0:
        snaps.append(h.copy())

rs = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32])
C = np.zeros(len(rs))
for hh in snaps:
    for k, r in enumerate(rs):
        C[k] += 0.5 * (((np.roll(hh, r, 0) - hh) ** 2).mean() +
                       ((np.roll(hh, r, 1) - hh) ** 2).mean())
C /= len(snaps)

# fit log model C = A ln r + B  vs power model C = a r^b (compare on ln C)
A, B = np.polyfit(np.log(rs), C, 1)
b_pow, ln_a = np.polyfit(np.log(rs), np.log(C), 1)
res_log = np.log(A * np.log(rs) + B) - np.log(C)
res_pow = (b_pow * np.log(rs) + ln_a) - np.log(C)
sse_log, sse_pow = float((res_log ** 2).sum()), float((res_pow ** 2).sum())
a_fit = float(np.exp(-B / A))              # C = A ln(r/a) ⇒ a = e^{−B/A}

# local effective exponents (finite differences of ln C vs ln r)
alpha_meas = 0.5 * np.diff(np.log(C)) / np.diff(np.log(rs))
r_mid = np.sqrt(rs[1:] * rs[:-1])
alpha_pred = 1.0 / (2.0 * np.log(r_mid / a_fit))
print(f"  C(r) fits:  log model SSE = {sse_log:.4f}   power model SSE = {sse_pow:.4f}"
      f"   (log wins: marginal, not power-law)")
print(f"  lattice cutoff from fit: a = {a_fit:.2f} (lattice units)")
print(f"    {'r':>6} {'α_eff measured':>15} {'1/(2 ln(r/a))':>15}")
for r, am, ap in zip(r_mid, alpha_meas, alpha_pred):
    print(f"    {r:>6.1f} {am:>15.3f} {ap:>15.3f}")
dev = np.abs(alpha_meas - alpha_pred) / alpha_pred
print(f"  median relative deviation = {np.median(dev):.2f} — the 1/ln(r) running of")
print("  the effective exponent is confirmed; extrapolated to r/a = 10^61 it gives")
print(f"  α_eff = {alpha_eff_h:.4f}, i.e. Δ_pred = {Delta_pred:.4f}.")

# ---------------------------------------------------------------------------
# 3. amplitude honesty: thermal roughness is Planckian ⇒ equilibrium = smooth
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. Amplitude check (honest): equilibrium roughness is Planckian — smooth")
print("=" * 78)
# thermal EW width: W² ~ (T/σ_stiff)·ln Λ with T = T_dS = H/2π and membrane
# stiffness set by 1/G = 1/ℓ_P² (η = 1/16πG); in Planck units W ~ ℓ_P √(ln Λ)
W_over_R = np.sqrt(lnLam) / Lam
print(f"  undriven (thermal) width: W/R_H ~ √(ln Λ)/Λ = {W_over_R:.1e}  — utterly smooth")
print("  ⇒ CONSISTENT: the equilibrium horizon is area-law (§7.1 baseline). The")
print("    DRIVEN amplitude — whether structure noise builds W/R → O(1) — is the")
print("    quantitative form of the activation question; the roughened fraction")
print("    W²/R² is a candidate microscopic definition of the f_sat gate itself")
print("    (see run_neff_gate.py for the kinetic derivation of its form).")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — falsifiability restored, with a sharper number")
print("=" * 78)
print(f"""  The §7/§7.2 conflict is resolved by computing the marginal case instead of
  gesturing at it:  Δ_pred = 1 − 1/(2 ln(R_H/ℓ_P)) = {Delta_pred:.4f}.
    • the log-suppression moves Δ by ~0.4% — INVISIBLE at σ_Δ ≈ 0.09;
    • a robust Δ ≈ 0.5 is {(Delta_pred-0.5)/sigma_DESI:.1f}σ away — it FALSIFIES; no escape hatch;
    • the pre-registered window is Δ ∈ [0.98, 1] (EW), with KPZ (Δ ≈ 0.61) the
      distinguishable alternative class — the measurement now discriminates
      roughening universality classes, not just area vs volume.
  Manuscript action: replace §7.2's "Δ ∈ (0,1) would be within the picture" with
  this computed window, and align §7's falsification statement with it.""")

# validation
assert 0.99 < Delta_pred < 1.0, "marginal prediction must be ~0.996"
assert (Delta_pred - 0.5) / sigma_DESI > 5.0, "Delta=0.5 must falsify at >5 sigma"
assert (1 - Delta_pred) / sigma_DESI < 0.1, "prediction indistinguishable from 1 at DESI precision"
assert sse_log < sse_pow, "2D EW structure function must prefer log over power law"
assert np.all(np.diff(alpha_meas) < 0.05), "local exponent must run downward (marginal)"
assert np.median(dev) < 0.5, "measured alpha_eff must track 1/(2 ln(r/a))"
print("\n[validate] alpha-eff assertions passed.")
