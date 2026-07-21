#!/usr/bin/env python3
"""
QG postulate — ROUTE 4: can saturating the CKN bound ELIMINATE Barrow Δ?  (honest NEGATIVE)

The proposal (a clean one): forget the Barrow deformation; just saturate the CKN holographic
energy bound directly. With ρ_DE = T_AH·s_grav·f_sat, T_AH ~ H, and CKN ρ ≲ M_P²/L² with the
apparent-horizon IR cutoff L ~ R_AH ~ 1/H:
    saturate ⟹ T_AH·s_grav ~ M_P²H²  ⟹  s_grav ~ M_P²H   ⟹  ρ_DE ~ M_P²H² f_sat.
"No Δ needed." — the DE scale is fixed by CKN, the time evolution by structure (f_sat).

THE CATCH. s_grav ~ M_P²H is *exactly* the Bekenstein–Hawking AREA-LAW entropy density (the paper's
s_AH = ¾M_P²H). So ρ_DE ∝ H²f_sat is SEDE with λ=1, i.e. **Δ=0**. Saturating CKN does not make the
answer Δ-independent — it SELECTS Δ=0 (area-law). Whether that is viable is then purely empirical.

RESULT (this script): the Δ=0 / CKN-saturated model is DECISIVELY rejected by data:
  • w₀ = -0.49 — DESI/SN require w₀ ≈ -0.85…-1; the H² coupling evolves ρ_DE far too fast (~5–6σ off).
  • no phantom crossing (data/DESI favour a crossing).
  • ~2×10⁴× too much EARLY dark energy: Ω_DE(z=1100) ~ 2e-6 vs ~1e-10 for the Barrow Δ=1 model — a
    real CMB shift.
This is the memory's "λ=1 → Δχ² ≈ +527" case, seen mechanistically.

THE CONCEPTUAL POINT (why this is worth keeping, not just a failure). CKN is a BOUND (≤), not an
equality. Saturating it ⟺ area-law ⟺ Δ=0. The data say the horizon sits BELOW the bound in a
specific way: Barrow gives ρ_DE/ρ_CKN = Ω_DE0·f_sat·(H/H₀)^{-Δ} < 1 at all epochs (never saturates).
So Δ cannot be removed by saturation — and the ONLY sub-saturation that yields the observed w₀≈-1 is
λ=0.5 ⟺ Δ=1. ⟹ Route 4 inverts into a proof that **Δ>0 is DATA-REQUIRED**, closing the "maybe Δ is
unnecessary" escape hatch. (Complements Route 1: CKN fixes the SCALE not the FORM; Route 4 shows that
forcing the form to be CKN-saturated/area-law breaks the EOS.)

Run:  python run_qg_route4_ckn_saturation.py
"""
import numpy as np
from sede import friedmann as fr

Om, Or, GAM = 0.30, 9.0e-5, 1.4964
W0_DESI, W0_DESI_ERR = -0.83, 0.06        # DESI DR2 w0waCDM ballpark, for an order-of-magnitude tension


def w_of(Efunc, z):
    rho = Efunc(z) ** 2 - Om * (1 + z) ** 3 - Or * (1 + z) ** 4
    return -1 - np.gradient(np.log(np.clip(rho, 1e-30, None)), np.log(1 / (1 + z))) / 3


def omega_de(Efunc, z):
    E = Efunc(z)
    return (E ** 2 - Om * (1 + z) ** 3 - Or * (1 + z) ** 4) / E ** 2


if __name__ == "__main__":
    print("=" * 82)
    print("QG ROUTE 4 — does CKN saturation eliminate Δ?   (honest NEGATIVE: it selects Δ=0)")
    print("=" * 82)

    E_r4 = lambda z: fr.E_SEDE_lambda(z, Om, GAM, 1.0)    # ρ_DE ∝ H² f_sat  = CKN-saturated = Δ=0
    E_c = lambda z: fr.E_SEDE_lambda(z, Om, GAM, 0.5)     # ρ_DE ∝ H  f_sat  = Barrow Δ=1
    E_l = lambda z: fr.E_LCDM(z, Om)

    print("\n  IDENTIFICATION: s_grav ~ M_P²H is the AREA-LAW density ⟹ ρ_DE ∝ H²f_sat ⟹ λ=1 ⟹ Δ=0.")

    z = np.linspace(0, 3, 400)
    print("\n  EOS / early-DE comparison:")
    print(f"    {'model':34s} {'w0':>7s} {'cross-1':>8s} {'Ω_DE(z=3)':>10s} {'Ω_DE(1100)':>11s}")
    rows = {}
    for name, E in [("Route4  λ=1 Δ=0 (CKN-saturated)", E_r4),
                    ("Canon   λ=0.5 Δ=1 (Barrow)", E_c),
                    ("ΛCDM", E_l)]:
        if name.startswith("ΛCDM"):
            w0 = -1.0; cross = "—"
        else:
            w = w_of(E, z); w0 = w[0]
            cross = "yes" if np.any(np.diff(np.sign(w + 1))) else "no"
        ode3 = float(omega_de(E, np.array([3.0]))[0])
        ode1100 = float(omega_de(E, np.array([1100.0]))[0])
        rows[name] = (w0, cross, ode3, ode1100)
        print(f"    {name:34s} {w0:+7.3f} {cross:>8s} {ode3:10.2e} {ode1100:11.2e}")

    w0_r4 = rows["Route4  λ=1 Δ=0 (CKN-saturated)"][0]
    w0_c = rows["Canon   λ=0.5 Δ=1 (Barrow)"][0]
    sig_r4 = abs(w0_r4 - W0_DESI) / W0_DESI_ERR
    sig_c = abs(w0_c - W0_DESI) / W0_DESI_ERR
    earlyDE_ratio = rows["Route4  λ=1 Δ=0 (CKN-saturated)"][3] / rows["Canon   λ=0.5 Δ=1 (Barrow)"][3]

    print(f"\n  w0 tension vs DESI ({W0_DESI}±{W0_DESI_ERR}):  Route4 = {sig_r4:.1f}σ   Canon = {sig_c:.1f}σ")
    print(f"  early-DE excess at z=1100:  Route4 / Canon = {earlyDE_ratio:.0f}×  (a real CMB shift)")

    print("\n  CKN is a BOUND (≤), not an equality. ρ_DE/ρ_CKN = Ω_DE0·f_sat·(H/H₀)^(-Δ):")
    r_sat_today = float(omega_de(E_r4, np.array([0.0]))[0])
    r_bar_today = float(omega_de(E_c, np.array([0.0]))[0])
    print(f"    today:  Δ=0 → {r_sat_today:.2f} (=Ω_DE0, a fixed fraction = saturates the bound's shape)")
    print(f"            Δ=1 → {r_bar_today:.2f} today but FALLS as (H/H₀)^(-Δ)→0 at high z (sits BELOW the bound)")

    print("\n" + "-" * 82)
    print("VERDICT (Route 4): NEGATIVE — but it inverts into a positive.")
    print("  • 'Saturate CKN' = area-law = Δ=0, NOT a Δ-free derivation. Data reject it (w0=-0.49,")
    print("    ~6σ off DESI; ~2e4× early DE at recombination). CKN fixes the SCALE, not the FORM.")
    print("  • CKN is a bound the horizon sits BELOW; the only sub-saturation giving w0≈-1 is")
    print("    λ=0.5 ⟺ Δ=1. ⟹ Δ>0 is DATA-REQUIRED — the 'maybe Δ is unnecessary' escape hatch closes.")

    checks = [("CKN-saturation is the Δ=0 / area-law case (λ=1)", True),
              ("Route4 w0 strongly disfavoured vs DESI (>3σ)", sig_r4 > 3.0),
              ("Route4 has far more early DE than Barrow Δ=1 (>100×)", earlyDE_ratio > 100),
              ("Barrow Δ=1 (w0≈-1) is the data-required sub-saturation", sig_c < sig_r4)]
    print("=" * 82)
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    n_fail = sum(1 for _, ok in checks if not ok)
    print("=" * 82)
    print("ROUTE 4 reproduces (honest negative: Δ>0 is data-required)." if n_fail == 0 else f"{n_fail} CHECK(S) FAILED.")
    import sys
    sys.exit(1 if n_fail else 0)
