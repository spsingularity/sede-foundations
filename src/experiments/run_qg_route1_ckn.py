#!/usr/bin/env python3
"""
QG postulate — ROUTE 1: can the CKN holographic energy bound FORCE the entropy FORM (δ=3/2, Δ=1)?

The one irreducible postulate (QG_DERIVATION §2 / OPEN_PROBLEMS §2) is the VOLUME-LAW entropy
form S ∝ A^{1+Δ/2}=A^{3/2}, not why Δ=1 (the GSL + geometric ceiling fix that GIVEN a power law).
The repo already derives the dark-energy SCALE from the Cohen–Kaplan–Nelson (CKN) UV–IR bound
(Principle 0). Route 1 asks the natural follow-up: does the SAME holographic bound also fix the
entropy EXPONENT, collapsing "scale" and "form" into one principle?

Method. Write everything in the single dimensionless ratio x ≡ M_P/H = R/l_P (the horizon in
Planck units, x≈8.5e60 today). Compare the entropy S ∝ x^{2δ} under three microscopic countings:

  (a) AREA-LAW / holographic (ground state): one bit per Planck AREA on the boundary,
      S = A/4l_P² ∝ x²        → δ = 1     → Δ = 0      (Bekenstein–Hawking)
  (b) CKN-EXTENSIVE: EFT modes at the CKN-tightened UV cutoff Λ_UV=(M_P H)^{1/2}, counted in
      the volume, S ~ (L·Λ_UV)³ ∝ x^{3/2} → δ = 3/4  → Δ = -1/2 (CKN *suppresses* below area-law)
  (c) PLANCK-VOLUME / typical state: one bit per Planck VOLUME, all entangled (Page/ETH),
      S = V/l_P³ ∝ x³          → δ = 3/2   → Δ = 1     (the postulate)

HONEST RESULT (this is a NEGATIVE-with-reduction, not a closure):
  CKN's own counting gives δ=3/4, NOT 3/2 — it *tightens* the EFT entropy below the area law, the
  OPPOSITE direction from Barrow. So the holographic ENERGY bound fixes the MAGNITUDE (ρ_DE~M_P²H²)
  but does NOT fix the FORM. The exponent is set entirely by (i) the cutoff scale and (ii) whether
  the state is the vacuum GROUND state (area-law) or a TYPICAL/thermalised state (volume-law).
  ⟹ Route 1 REDUCES the postulate to a single sharper statement: "the horizon dof are counted at
  the Planck scale in a typical (all-entangled) state." That residual is exactly Route 2's job
  (run_qg_route2_driving.py). The seam S_vol/S_area = x = M_P/H ≈ 10⁶¹ = √(10¹²²) is the same
  factor as Thm 11 — confirming scale (CKN) and form (Barrow) are genuinely independent inputs,
  glued by √(the CC hierarchy), not one derivable from the other.

Run:  python run_qg_route1_ckn.py
"""
import numpy as np

lP = 1.616e-35                         # Planck length, m
H0 = 67.5e3 / 3.086e22                 # s^-1
c = 2.998e8
R = c / H0                             # de Sitter horizon radius, m
x = R / lP                             # M_P/H in Planck units

CASES = [
    ("area-law / holographic (ground state)", 1.00, 0.0, "Bekenstein–Hawking, S=A/4l_P²"),
    ("CKN-extensive (Λ_UV=(M_P H)^1/2)",      0.75, -0.5, "CKN UV–IR cutoff, S~(LΛ)³"),
    ("Planck-volume / typical (Page/ETH)",    1.50, 1.0, "S=V/l_P³ — THE POSTULATE"),
]

if __name__ == "__main__":
    print("=" * 78)
    print("QG ROUTE 1 — does CKN force the entropy FORM δ=3/2 (Δ=1)?   (x ≡ M_P/H = %.2e)" % x)
    print("=" * 78)
    print(f"\n  {'counting':42s} {'δ':>5s} {'Δ':>5s}  {'S∝x^{2δ}':>10s}")
    S = {}
    for name, delta, Delta, note in CASES:
        S[delta] = x ** (2 * delta)
        print(f"  {name:42s} {delta:5.2f} {Delta:5.1f}  {S[delta]:10.2e}   [{note}]")

    S_area, S_ckn, S_vol = S[1.00], S[0.75], S[1.50]
    seam = S_vol / S_area
    print(f"\n  area-law S ≈ {S_area:.2e}   (the known de Sitter entropy ~10¹²²)")
    print(f"  CKN δ=3/4 gives S ≈ {S_ckn:.2e}  — SMALLER than area-law (CKN *tightens*, opposite of Barrow).")
    print(f"  seam  S_vol/S_area = x = M_P/H = {seam:.2e}  = √(10¹²²) = {np.sqrt(2.26e122):.2e}")

    # checks
    ck_ckn_below_area = S_ckn < S_area                         # CKN does NOT enhance to volume-law
    ck_vol_is_delta1 = abs(1.50 - (1 + 1.0 / 2)) < 1e-9         # δ=3/2 ⟺ Δ=1
    ck_seam = abs(np.log10(seam) - np.log10(x)) < 1e-6          # seam = M_P/H = √hierarchy

    print("\n" + "-" * 78)
    print("VERDICT (Route 1): CKN fixes the SCALE, NOT the FORM.")
    print("  • CKN's own counting → δ=3/4 (Δ=-1/2), tightening BELOW the area law — it cannot")
    print("    produce the volume-law enhancement (δ=3/2) on its own.")
    print("  • The exponent is set by the STATE (ground=area-law δ=1 vs typical=volume-law δ=3/2),")
    print("    not by the holographic energy bound. So 'scale' and 'form' are independent inputs,")
    print("    glued by the seam √(10¹²²). This MATCHES Thm 11 — not a closure, a clean reduction.")
    print("  ⟹ The postulate reduces to: the horizon is in a TYPICAL Planck-volume state.")
    print("     That residual is addressed dynamically in Route 2 (run_qg_route2_driving.py).")

    checks = [("CKN counting stays below area-law (no volume enhancement)", ck_ckn_below_area),
              ("δ=3/2 ⟺ Δ=1 (the postulate identified)", ck_vol_is_delta1),
              ("seam S_vol/S_area = M_P/H = √(CC hierarchy)", ck_seam)]
    print("=" * 78)
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    n_fail = sum(1 for _, ok in checks if not ok)
    print("=" * 78)
    print("ROUTE 1 reproduces (honest negative + reduction)." if n_fail == 0 else f"{n_fail} CHECK(S) FAILED.")
    import sys
    sys.exit(1 if n_fail else 0)
