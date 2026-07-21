"""
run_gr_height.py — R4 at derivation level: the linearised-GR ingredients of the
horizon height equation, computed exactly where possible.

THE GAP. run_dns_roughening.py assembled ∂_t h = D∇²h − κh + η_dep from
order-of-magnitude membrane ingredients (κ ≈ H). This script replaces the
order-of-magnitude steps with exact GR results and recomputes the dissipation
band ε(z), leaving one clearly-scoped piece (the full ℓ-decomposed FOTH
perturbation) as the remaining derivation.

FOUR RESULTS:
  (1) EXACT SINK RATE (symbolic). The Cai–Kim surface gravity of the flat-FRW
      apparent horizon, derived from its definition
      κ = (1/2√−h) ∂_a(√−h h^{ab} ∂_b R̃) at R̃ = 1/H, is
          κ = −H·(1 − 3w_eff)/4,
      verified with sympy from the 2-metric. Limits: matter era |κ| = H/4,
      de Sitter |κ| = H, radiation κ = 0. The sink is SLOWER than the H used in
      the prototype — Fork B is strengthened.
  (2) EXACT DEPOSITION RESPONSE (symbolic). Perturbing the Misner–Sharp
      apparent-horizon condition 2GM(R_A) = R_A with a mass deposit δM gives
          δR_A/R_A = −(1/2)·δM/M  —
      the source term of the height equation, with an exact coefficient (the
      horizon responds linearly to deposited mass-energy; sign: overdensity ⇒
      apparent horizon moves inward).
  (3) THE SINK ACTS ONLY AT HORIZON-SCALE MODES (estimate + exact κ). With the
      FOTH transport coefficients (η = ζ = +1/16πG for dynamical horizons,
      Gourgoulhon–Jaramillo 2006 — positive, unlike the teleological event
      horizon) and membrane surface density Σ ∼ κ/8πG, the height diffusivity is
      D ∼ η/Σ = 1/(2κ), so the sink competes with redistribution only for
      k ≲ k_c = √2·κ: k_c·R_A = √2(1−3w)/4 ≲ 1. Sub-horizon modes are
      redistribution-dominated — "dissipation only at the Hubble boundary"
      (Hwa–Kardar condition) realised as boundary-in-mode-space, now with the
      exact κ.
  (4) THE UPDATED ε(z) BAND. ε(z) = |κ(z)|·t_dyn = (1−3w_eff(z))/4 · Δ_vir^{−1/2}
      on the actual SEDE background: ε ∈ [~0.02, ~0.06] over the structure era —
      SMALLER than the 0.075 prototype value, hence a LONGER truncated-SOC
      window, s_c = 2/ε² ∈ [~10³, ~10⁴].

REMAINING (scoped, not hidden): the transverse diffusivity D is an estimate
(η/Σ); the complete derivation is the ℓ-decomposed linear perturbation of a
future outer trapping horizon in FRW, which would fix D's coefficient. κ (the
sink) and the source coefficient are now exact; the CLASS conclusion (nearly
conservative, ε ≪ 1) no longer rests on order-of-magnitude input.
"""
import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# 1. exact Cai–Kim surface gravity from the definition (sympy)
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. EXACT sink rate: κ = −H(1 − 3w)/4 from the Cai–Kim definition (sympy)")
print("=" * 78)
t, r, w = sp.symbols("t r w", positive=True)
a = sp.Function("a", positive=True)(t)
H = sp.diff(a, t) / a
Rt = a * r                                            # areal radius R̃ = a(t)·r
# 2-metric h_ab = diag(−1, a²): √−h = a, h^tt = −1, h^rr = 1/a²
sqh = a
term_t = sp.diff(sqh * (-1) * sp.diff(Rt, t), t)      # ∂_t(√−h h^tt ∂_t R̃)
term_r = sp.diff(sqh * (1 / a**2) * sp.diff(Rt, r), r)
kappa = sp.simplify((term_t + term_r) / (2 * sqh))
# evaluate on the apparent horizon r_A = 1/ȧ  (R̃ = 1/H)
kappa_AH = sp.simplify(kappa.subs(r, 1 / sp.diff(a, t)))
# substitute the Friedmann dynamics: ä = −(1/2)(1+3w)H²a  (flat, EoS w)
addot = -(sp.Rational(1, 2)) * (1 + 3 * w) * H**2 * a
kappa_w = sp.simplify(kappa_AH.subs(sp.diff(a, t, 2), addot))
kappa_pred = -H * (1 - 3 * w) / 4
diff_expr = sp.simplify(kappa_w - kappa_pred)
print(f"  κ(AH) from definition   : {sp.simplify(kappa_AH)}")
print(f"  with ä = −(1+3w)H²a/2   : κ = {kappa_w}")
print(f"  claimed closed form     : κ = −H(1−3w)/4;  difference = {diff_expr}")
for name, wv, expect in [("radiation", sp.Rational(1, 3), 0),
                         ("matter", 0, sp.Rational(1, 4)),
                         ("de Sitter", -1, 1)]:
    val = sp.simplify(-kappa_w.subs(w, wv) / H)
    print(f"    {name:>10}: |κ|/H = {val}   (expected {expect})")
print("  ⇒ the prototype's κ ≈ H was an OVERestimate: the exact sink is H(1−3w)/4,")
print("    between H/4 (matter era) and H (de Sitter) — the horizon is MORE nearly")
print("    conservative than assumed.")

# ---------------------------------------------------------------------------
# 2. exact deposition response from the Misner–Sharp condition (sympy)
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. EXACT source: δR_A/R_A = −(1/2)·δM/M from the Misner–Sharp condition")
print("=" * 78)
G, R, dR, dM, rho = sp.symbols("G R deltaR deltaM rho", positive=True)
M_of_R = sp.Rational(4, 3) * sp.pi * rho * R**3       # background mass inside R
# AH condition: 2G·M(R_A) = R_A  (⇒ R_A = 1/H, since 8πGρ/3 = H²)
cond0 = sp.Eq(2 * G * M_of_R, R)
rho_sol = sp.solve(cond0, rho)[0]                     # ρ on the AH
# perturb: 2G[M(R+δR) + δM] = R + δR, linearise in δR (Newton form)
F = (2 * G * (M_of_R.subs(R, R + dR) + dM) - (R + dR)).subs(rho, rho_sol)
dR_sol = sp.simplify(-F.subs(dR, 0) / sp.diff(F, dR).subs(dR, 0))
M_AH = M_of_R.subs(rho, rho_sol)                      # = R/2G
ratio = sp.simplify((dR_sol / R) / (dM / M_AH))
print(f"  linearised solution: δR_A = {sp.simplify(dR_sol)}")
print(f"  ⇒ (δR_A/R_A)/(δM/M) = {ratio}   (exact: −1/2)")
print("  The horizon responds LINEARLY to deposited mass-energy with an exact O(1)")
print("  coefficient — the deposition source of the height equation is derived, not")
print("  modelled; sign: an overdensity pulls the apparent horizon inward.")

# ---------------------------------------------------------------------------
# 3. the sink is confined to horizon-scale modes (FOTH transport + exact κ)
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. Mode-space locality of the sink: k_c·R_A = √2·|κ|/H ≲ 1")
print("=" * 78)
print("  FOTH transport (Gourgoulhon–Jaramillo 2006): η = ζ = +1/16πG (positive —")
print("  diffusive, unlike the teleological event horizon's ζ = −1/16πG).")
print("  With membrane surface density Σ ∼ κ/8πG the height diffusivity is")
print("  D ∼ η/Σ = 1/(2κ), so sink beats redistribution only for k < k_c = √2·κ:")
for name, wv in [("matter era", 0.0), ("today (w_eff≈−0.7)", -0.7), ("de Sitter", -1.0)]:
    kH = (1 - 3 * wv) / 4
    print(f"    {name:>20}: |κ|/H = {kH:.2f}   k_c·R_A = {np.sqrt(2)*kH:.2f}")
print("  ⇒ every sub-horizon mode is redistribution-dominated; the only mode the")
print("    sink controls is the horizon-scale (Hubble) mode — the Hwa–Kardar")
print("    'dissipation at the boundary only' condition realised in mode space,")
print("    now with the exact κ. [D's O(1) coefficient remains the scoped residual.]")

# ---------------------------------------------------------------------------
# 4. the updated ε(z) band on the actual SEDE background
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("4. ε(z) = |κ(z)|·t_dyn on the SEDE background — the band, and the new cutoff")
print("=" * 78)
from sede.friedmann import E_SEDE_lambda
OM, GAM, LAM = 0.30, 1.4964, 0.5
zg = np.linspace(0.0, 5.0, 400)
E = E_SEDE_lambda(zg, OM, GAM, LAM)
dlnE_dz = np.gradient(np.log(E), zg)
w_eff = -1.0 + (2.0 / 3.0) * (1 + zg) * dlnE_dz       # total effective EoS
kap_H = (1 - 3 * w_eff) / 4
DVIR = 178.0
eps_z = kap_H * DVIR ** -0.5
print(f"    {'z':>5} {'w_eff':>8} {'|κ|/H':>7} {'ε(z)':>8} {'s_c=2/ε²':>9}")
for zz in [0.0, 0.5, 1.0, 2.0, 5.0]:
    k = np.argmin(np.abs(zg - zz))
    print(f"    {zz:>5.1f} {w_eff[k]:>8.2f} {kap_H[k]:>7.2f} {eps_z[k]:>8.3f} {2/eps_z[k]**2:>9.0f}")
band = (float(eps_z[(zg > 0) & (zg < 3)].min()), float(eps_z[(zg > 0) & (zg < 3)].max()))
print(f"  structure era (0 < z < 3): ε ∈ [{band[0]:.3f}, {band[1]:.3f}]  ⇒  s_c ∈ "
      f"[{2/band[1]**2:.0f}, {2/band[0]**2:.0f}]")
print("  ⇒ the exact κ LOWERS ε below the 0.075 prototype value throughout the")
print("    structure era: the truncated-SOC window is LONGER (3–4 decades), and the")
print("    near-conservation of the roughening is now anchored in exact GR")
print("    quantities (κ, source) rather than order-of-magnitude estimates.")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — R4 upgraded from order-of-magnitude to derivation level")
print("=" * 78)
print(f"""  Derived exactly (sympy): the sink rate κ = −H(1−3w)/4 from the Cai–Kim
  definition (matter H/4, dS H, radiation 0), and the deposition source
  δR_A/R_A = −(1/2)δM/M from the Misner–Sharp condition. Cited (established):
  the FOTH transport coefficients η = ζ = +1/16πG. Estimated (scoped residual):
  the height diffusivity D ∼ 1/(2κ), pending the ℓ-decomposed FOTH perturbation.
  Consequences: ε(z) = (1−3w_eff)/4·Δ_vir^(−1/2) ∈ [{band[0]:.2f}, {band[1]:.2f}] over the
  structure era — smaller than the 0.075 prototype ⇒ a longer τ = 3/2 window
  (s_c up to ~{2/band[0]**2:.0f}); the sink acts only at horizon-scale modes (k_c·R_A ≲ 1),
  which IS the Hwa–Kardar boundary-dissipation condition. The dns-roughening
  ledger row can move from ARGUED (prototype) to DERIVED (κ, source exact;
  D coefficient scoped) — the class conclusion no longer rests on κ ≈ H.""")

# validation
assert diff_expr == 0, "Cai-Kim kappa must equal -H(1-3w)/4 symbolically"
assert sp.simplify(-kappa_w.subs(w, 0) / H - sp.Rational(1, 4)) == 0, "matter limit H/4"
assert sp.simplify(-kappa_w.subs(w, -1) / H - 1) == 0, "de Sitter limit H"
assert sp.simplify(kappa_w.subs(w, sp.Rational(1, 3))) == 0, "radiation limit 0"
assert sp.simplify(ratio + sp.Rational(1, 2)) == 0, "MS response must be -1/2 exactly"
assert 0.01 < band[0] < band[1] < 0.08, "epsilon band must be small (near-conservative)"
assert float(np.sqrt(2) * kap_H[(zg > 0) & (zg < 3)].max()) < 1.3, \
    "sink confined to horizon-scale modes over the structure era"
print("\n[validate] gr-height assertions passed.")
