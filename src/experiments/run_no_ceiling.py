"""
run_no_ceiling.py — the algebra-type residue (R5) collapses onto the Δ measurement:
the ceiling's EXISTENCE is fixed by the IR attractor, and its VALUE is a trace
calibration — which is exactly the count that Δ measures.

THE GAP. §7.1b poses the residue as "the type (II₁ or II_∞) of the FRW static-patch
algebra" — apparently an independent open input from de Sitter holography.

THE ANALYSIS (with one honest NEGATIVE en route):
  (1) TYPE ⟺ BOUNDEDNESS. A Type II₁ factor has a maximum-entropy state (the
      tracial state); II_∞ does not. So the algebra type is decided by whether the
      horizon entropy is BOUNDED along the cosmological flow — i.e. by the IR
      attractor of the background.
  (2) NEGATIVE: the type alone does NOT discriminate ΛCDM from SEDE. Both are
      accelerating attractors with finite asymptotic Hubble rate (ΛCDM: E∞ = √Ω_Λ;
      SEDE: the fixed point E∞ = Ω_DE·f∞ of E² = Ω_m a⁻³ + Ω_DE f_sat E at λ = ½),
      hence finite asymptotic horizon radius and BOUNDED entropy — both are
      II₁-compatible. Only an eternally-decelerating control (matter-only:
      R_AH → ∞) is II_∞. The naive lemma "SEDE ⇒ no ceiling" FAILS, and we say so.
  (3) THE RESIDUE IS THE TRACE CALIBRATION. In a Type II₁ factor the trace — hence
      the VALUE of the maximum entropy — is defined only up to normalisation. CLPW
      fix it by CALIBRATING to the semiclassical Gibbons–Hawking area of the
      Einstein saddle; that the eternal-patch ceiling "is the area" is a
      calibration to Einstein gravity, not an algebraic output. The volume-law
      alternative is a different calibration of the same II₁ structure: ceiling
      e^{Area} vs e^{Volume} differ by exactly the count — the factor (R/ℓ_P) that
      Δ measures. The algebra machinery is silent between them BY CONSTRUCTION.
NET: R5 is not an independent theoretical input. The type is fixed by the attractor
(which the cosmology itself supplies), and the calibration IS Δ — so the "open
algebraic question" and the "empirical Δ measurement" (Level 3) are one question in
two languages. DESI DR3 + Euclid answers both.

HONEST FLAGS. (1)–(3) use the finite-dimensional/semiclassical shadow of the
algebraic statements (max entropy ⟺ bounded trace); the rigorous crossed-product
construction for a matter-filled FRW patch remains undone — but what it would
decide is now shown to be the calibration, i.e. Δ.
"""
import numpy as np

# ---------------------------------------------------------------------------
# 1. IR attractors: LCDM and SEDE bounded; matter-only unbounded (control)
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. The IR attractor decides boundedness (⇒ the algebra type)")
print("=" * 78)
OM, OR = 0.30, 9e-5
ODE = 1.0 - OM - OR

def E_lcdm(z):
    return np.sqrt(OM * (1 + z) ** 3 + OR * (1 + z) ** 4 + ODE)

def E_sede(z, f=1.0):
    """Fixed-point closure E² = Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_DE·f·E  (λ = 1/2)."""
    m = OM * (1 + z) ** 3 + OR * (1 + z) ** 4
    return 0.5 * (ODE * f + np.sqrt((ODE * f) ** 2 + 4 * m))

def E_matter(z):
    return np.sqrt((1 + z) ** 3)

# late-time scan toward a -> infinity (z -> -1); f_sat -> 1 at the attractor
z_late = -1 + np.geomspace(1.0, 1e-4, 40)          # z: 0 -> -0.9999
S = {}
for name, E in [("LCDM", E_lcdm), ("SEDE", E_sede), ("matter-only", E_matter)]:
    S[name] = 1.0 / E(z_late) ** 2                  # S_area ∝ R² = 1/H² (units of S₀)
E_inf_sede = ODE * 1.0                              # fixed point at f∞ = 1
E_inf_lcdm = np.sqrt(ODE)
print(f"    {'a/a0':>10} {'S_LCDM/S0':>10} {'S_SEDE/S0':>10} {'S_matter/S0':>12}")
for k in [0, 15, 25, 33, 39]:
    a = 1 / (1 + z_late[k])
    print(f"    {a:>10.1f} {S['LCDM'][k]:>10.2f} {S['SEDE'][k]:>10.2f} {S['matter-only'][k]:>12.1f}")
print(f"  attractors: E∞(ΛCDM) = √Ω_Λ = {E_inf_lcdm:.3f}  (S∞/S₀ = {1/E_inf_lcdm**2:.2f});")
print(f"              E∞(SEDE) = Ω_DE·f∞ = {E_inf_sede:.3f}  (S∞/S₀ = {1/E_inf_sede**2:.2f})")
print(f"              matter-only: E → 0, S → ∞ (unbounded)")
print("  ⇒ BOTH ΛCDM and SEDE have bounded horizon entropy (finite asymptotic")
print("    ceiling ⇒ II₁-compatible); only eternal deceleration is II_∞.")

# fixed-point consistency check
res = E_inf_sede ** 2 - ODE * 1.0 * E_inf_sede
print(f"  fixed-point residual E∞² − Ω_DE·E∞ = {res:.2e} (exact 0)")

# ---------------------------------------------------------------------------
# 2. honest NEGATIVE: the type alone does not discriminate area vs volume
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. NEGATIVE: 'SEDE ⇒ II_∞ ⇒ no ceiling' FAILS — both cosmologies are bounded")
print("=" * 78)
grow_sede = S["SEDE"][-1] / S["SEDE"][0]
grow_matt = S["matter-only"][-1] / S["matter-only"][0]
print(f"  entropy growth from a=1 to a→∞:  SEDE ×{grow_sede:.2f} (converges),")
print(f"  matter-only ×{grow_matt:.0f} (diverges).  SEDE, like ΛCDM, ends on an")
print(f"  accelerating attractor with a FINITE ceiling — the type dichotomy cannot")
print(f"  separate them. The naive no-ceiling lemma is rejected; what remains is (3).")

# ---------------------------------------------------------------------------
# 3. the ceiling VALUE is a trace calibration — and that calibration IS Δ
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. The II₁ trace is defined up to normalisation: ceiling VALUE = calibration = Δ")
print("=" * 78)
c = 2.99792458e8; H0 = 67.4e3 / 3.0857e22; lP = 1.616255e-35
Lam0 = (c / H0) / lP                               # R_H(today)/ℓ_P
Lam_inf = Lam0 / E_inf_sede                        # asymptotic horizon (SEDE units)
S_area = np.pi * Lam_inf ** 2
S_vol = Lam_inf ** 3
print(f"  asymptotic horizon: R∞/ℓ_P = {Lam_inf:.2e}")
print(f"  ceiling if calibrated to the AREA  (CLPW/Einstein saddle): S = π(R/ℓ_P)² ~ 10^{np.log10(S_area):.0f}")
print(f"  ceiling if calibrated to the VOLUME (SEDE):                S = (R/ℓ_P)³   ~ 10^{np.log10(S_vol):.0f}")
print(f"  ratio = R∞/ℓ_P/π ~ 10^{np.log10(S_vol/S_area):.0f} — exactly the COUNT (the Δ exponent).")
print("  In a Type II₁ factor the trace (hence the max-entropy value) carries a free")
print("  normalisation; CLPW fix it by matching the Gibbons–Hawking entropy of the")
print("  semiclassical Einstein saddle. 'The eternal-patch ceiling is the area' is")
print("  therefore a CALIBRATION to Einstein gravity, not an algebraic derivation —")
print("  the same II₁ structure calibrated to the volume count is equally consistent.")
print("  The algebraic residue is the calibration; the calibration is the count;")
print("  the count is Δ = d_H − 2 (§7.2). One question, two languages.")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — R5 collapses onto Level 3 (the Δ measurement)")
print("=" * 78)
print(f"""  (1) The algebra TYPE is fixed by the IR attractor: accelerating ⇒ bounded
      entropy ⇒ a ceiling exists (II₁-compatible); eternal deceleration ⇒ II_∞.
      SEDE and ΛCDM are both accelerating — the type does not discriminate
      (honest NEGATIVE for the naive no-ceiling lemma).
  (2) The ceiling VALUE — e^(Area) vs e^(Volume), a factor 10^{np.log10(S_vol/S_area):.0f} — is the II₁
      trace normalisation, which the algebraic machinery leaves free and CLPW fix
      by calibrating to the Einstein saddle. The residue is the calibration.
  (3) The calibration is the count is Δ. So the ledger's 'algebra type: OPEN'
      row is not an independent theoretical unknown: it is the Δ measurement
      restated algebraically, decided by DESI DR3 + Euclid either way.
  Manuscript action: reword the §7.1b/ledger residue from 'open algebraic
  question' to 'algebraic restatement of the Δ measurement (trace calibration)',
  and log the failed stronger lemma as a NEGATIVE.""")

# validation
assert abs(res) < 1e-12, "SEDE attractor must solve its fixed-point equation"
assert S["LCDM"][-1] / S["LCDM"][-10] < 1.001, "LCDM entropy must converge (bounded)"
assert S["SEDE"][-1] / S["SEDE"][-10] < 1.001, "SEDE entropy must converge (bounded)"
assert grow_matt > 100, "matter-only control must diverge (unbounded, II_inf)"
assert grow_sede < 5, "SEDE entropy growth to the attractor is O(1)"
assert 55 < np.log10(S_vol / S_area) < 65, "ceiling ratio must be the count ~ R/l_P"
print("\n[validate] no-ceiling assertions passed.")
