#!/usr/bin/env python3
"""
QG postulate — ROUTE 3: the horizon as a roughening surface (Δ as a universality class).

Borrowed + merged from the sibling SEDE_V2 team's S1/KPZ work (theory/QG_POSTULATE_STRATEGIES.md,
run_kpz_roughening.py), reimplemented here and connected to our Routes 1 & 2.

THE REFRAME (their contribution). Stop counting static microstates; ask "why does the cosmic horizon
become space-filling?". A volume-law entropy in embedding dimension d gives S∝R^d, ordinary area
A∝R^{d-1}, hence S∝A^{1+1/(d-1)} ⟹ Δ = 2/(d-1). So Δ=1 is the d=3 (space-filling) ENDPOINT, not an
awkward intermediate — dissolving the "why this Δ" objection. [This is the Route-1 'volume-law'
exponent in geometric language; cf. run_qg_route1_ckn.py.]

THE DYNAMICS (their KPZ upgrade). Model the horizon as a stochastic growing 2-surface,
    ∂h/∂t = ν∇²h + (λ/2)(∇h)² + η ,
a self-affine surface whose graph has Hausdorff dimension d_H = 3 - α (α = roughness exponent), so
Barrow's deformation (its ORIGINAL fractal-area definition) is Δ = d_H - 2 = 1 - α. Then Δ is fixed
by the UNIVERSALITY CLASS of the roughening:
  • Edwards-Wilkinson (λ=0, linear/thermal): 2+1D is lower-critical, α→0 (log) ⟹ Δ→1  [TEXTBOOK].
  • KPZ (λ≠0, nonlinear lateral growth): relevant in 2+1D, α≈0.39 ⟹ Δ≈0.6        [TEXTBOOK].
A finite-N integration here reproduces only the TREND (stronger λ ⟹ larger α ⟹ smaller Δ); the
asymptotic values above are established universality theory, not measured from this small grid.

NOTE on the two Δ-relations: Δ=2/(d-1) uses the VOLUME embedding dimension d (=3 fixed); Δ=d_H-2
uses the SURFACE Hausdorff dimension d_H (∈[2,3]). They are DIFFERENT parametrizations that coincide
only at the space-filling endpoint d=d_H=3 ⟹ Δ=1 (the robust shared conclusion). We keep them
distinct on purpose.

THE MERGE (our contribution). The stochastic noise η that drives the roughening is, in SEDE, the
structure-sourced entropy deposition of Route 2 (run_qg_route2_driving.py) — and their strategy S4
("maximise entropy PRODUCTION / driven steady state") IS our Route 2 (independent arrival). The
EW class condition (λ=0, no nonlinear lateral-growth term) is exactly SEDE's founding postulate
F = ρ - T·s = 0 (purely entropic; Thm 9 Pillar 1): with zero free energy there is no lateral
'growth' drive, only thermal roughening ⟹ EW, not KPZ. So:

    structure-driving (Route 2) = the noise η   →   EW-class roughening   →   Δ = 1,
    and "why EW not KPZ" reduces to SEDE's already-stated F=0 postulate (data favour it, below).

This does NOT close the postulate (the QG roughening LAW, and why the cosmic horizon is thermal/rough,
stay open) — but it ties Δ=1 to an assumption SEDE ALREADY makes, with a falsifiable data preference.

Run:  python run_qg_route3_roughening.py
"""
import numpy as np

RNG = np.random.default_rng(3)
DELTA_DATA, DELTA_ERR = 0.93, 0.09          # cross-team Barrow Δ (MCMC); our forecast σ(Δ)≈0.09 consistent


def _lap(h):
    return np.roll(h, 1, 0) + np.roll(h, -1, 0) + np.roll(h, 1, 1) + np.roll(h, -1, 1) - 4 * h


def _grad2(h):
    gx = 0.5 * (np.roll(h, -1, 0) - np.roll(h, 1, 0))
    gy = 0.5 * (np.roll(h, -1, 1) - np.roll(h, 1, 1))
    return gx ** 2 + gy ** 2


def _evolve(N, nu, lam, D, dt, steps):
    h = np.zeros((N, N))
    amp = np.sqrt(2 * D * dt)
    for _ in range(steps):
        h += dt * (nu * _lap(h) + 0.5 * lam * _grad2(h)) + amp * RNG.standard_normal((N, N))
        h -= h.mean()
    return h


def _alpha(h):
    """Roughness exponent from the structure function C(r)=<(h(x+r)-h(x))²> ~ r^{2α}."""
    N = h.shape[0]
    rs = np.unique(np.round(np.geomspace(1, N // 4, 10)).astype(int))
    C = [((np.roll(h, r, 0) - h) ** 2 + (np.roll(h, r, 1) - h) ** 2).mean() / 2.0 for r in rs]
    return float(max(0.0, np.polyfit(np.log(rs), np.log(C), 1)[0] / 2.0))


def _run(label, N, nu, lam, D, dt, steps):
    a = _alpha(_evolve(N, nu, lam, D, dt, steps))
    Delta = 1.0 - a                                  # Δ = d_H - 2 = 1 - α
    print(f"    [{label:20s}] α={a:.3f}  d_H=3-α={3-a:.3f}  Δ=1-α={Delta:.3f}")
    return Delta


if __name__ == "__main__":
    print("=" * 82)
    print("QG ROUTE 3 — horizon roughening: Δ as a universality class (borrowed from SEDE_V2 S1/KPZ)")
    print("=" * 82)

    # (1) the geometric reframe: Δ = 2/(d-1) -> Δ=1 is the d=3 space-filling endpoint
    print("\n  (1) Volume-law dimension relation  Δ = 2/(d-1)  (embedding dim d):")
    for d in (3, 4, 6, 100):
        print(f"        d={d:>3d}: Δ={2/(d-1):.3f}" + ("   <- SEDE: space-filling in 3-space" if d == 3 else ""))
    endpoint_ok = abs(2 / (3 - 1) - 1.0) < 1e-9

    # (2) the dynamical universality class (finite-N TREND only; asymptotics are textbook)
    print("\n  (2) Roughening flow ∂h/∂t=ν∇²h+(λ/2)(∇h)²+η  (128², finite-N — TREND, not asymptote):")
    N, D, dt, steps = 128, 1.0, 0.10, 6000
    d_ew = _run("Edwards-Wilkinson λ=0", N, 1.0, 0.0, D, dt, steps)
    d_k2 = _run("KPZ λ=2", N, 1.0, 2.0, D, 0.04, steps)
    d_k4 = _run("KPZ λ=4", N, 1.0, 4.0, D, 0.02, steps)
    trend_ok = d_ew > d_k2 > d_k4                     # stronger nonlinearity lowers Δ
    print(f"    trend: EW {d_ew:.2f} > KPZ-λ2 {d_k2:.2f} > KPZ-λ4 {d_k4:.2f}  "
          f"(stronger KPZ nonlinearity -> larger α -> smaller Δ)")
    print("    asymptotic universality (textbook): EW 2+1D lower-critical α→0 ⟹ Δ→1 ; KPZ α≈0.39 ⟹ Δ≈0.6")

    # (3) which class do the data pick?
    s_ew = abs(1.0 - DELTA_DATA) / DELTA_ERR
    s_kpz = abs(0.6 - DELTA_DATA) / DELTA_ERR
    print(f"\n  (3) Data Δ={DELTA_DATA}±{DELTA_ERR}:  vs EW(Δ→1) = {s_ew:.1f}σ ;  vs KPZ(Δ≈0.6) = {s_kpz:.1f}σ")
    print(f"      -> data FAVOUR the EW/thermal (F=0) class over generic KPZ.")
    data_ok = s_ew < s_kpz

    print("\n" + "-" * 82)
    print("THE MERGE (Route 2 ⟷ Route 3 ⟷ SEDE's F=0 postulate):")
    print("  • The noise η driving the roughening = SEDE's structure-sourced entropy deposition")
    print("    (Route 2). Their strategy S4 (max entropy PRODUCTION / driven steady state) = our Route 2.")
    print("  • EW class (λ=0, no nonlinear lateral-growth term) ⟺ SEDE's founding F=ρ-Ts=0 (Thm 9 P1):")
    print("    zero free energy ⟹ no lateral growth drive ⟹ thermal (EW) roughening ⟹ Δ→1.")
    print("  ⟹ 'why volume-law Δ=1' reduces to 'why EW not KPZ' reduces to SEDE's ALREADY-STATED F=0,")
    print("     with the data preferring it (EW 0.8σ vs KPZ 3.7σ). NOT a closure — the QG roughening")
    print("     law and the horizon's thermal/rough STATE remain the one open postulate (OPEN_PROBLEMS §2).")

    checks = [("Δ=2/(d-1): d=3 is the space-filling Δ=1 endpoint", endpoint_ok),
              ("roughening trend: stronger KPZ nonlinearity lowers Δ", trend_ok),
              ("data favour EW (Δ→1) over KPZ (Δ≈0.6)", data_ok)]
    print("=" * 82)
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    n_fail = sum(1 for _, ok in checks if not ok)
    print("=" * 82)
    print("ROUTE 3 reproduces (borrowed + merged with Routes 1/2)." if n_fail == 0 else f"{n_fail} CHECK(S) FAILED.")
    import sys
    sys.exit(1 if n_fail else 0)
