"""
run_dns_roughening.py — deciding Fork A vs Fork B from horizon hydrodynamics:
the membrane equations put the horizon in the NEARLY-CONSERVATIVE class
automatically, with a computable dissipation parameter ε = H·t_dyn ≈ 0.075.

THE GAP. The ledger's one OPEN dynamical row: "conservation law — bulk-conserved
roughening current, boundary-only loss — decides Fork A vs B." The paper leaves
this as an unknown property the horizon may or may not have.

THE CLOSURE. The horizon already HAS first-principles hydrodynamics — we don't
need to postulate the roughening equation:
  (1) TRANSPORT (Damour–Navier–Stokes): the membrane momentum equation is a LOCAL
      CONSERVATION LAW with known coefficients; for the dynamical/apparent
      horizon (SEDE's choice) the bulk viscosity is POSITIVE, ζ = +1/16πG
      (Gourgoulhon–Jaramillo 2005), unlike the teleological event horizon's
      ζ = −1/16πG — so deformation REDISTRIBUTES diffusively (ν∇²h = −∇·J_h,
      locally conserving), it does not anti-diffuse.
  (2) SINK (linearised Raychaudhuri): deformations of the apparent horizon relax
      at the surface-gravity rate κ ≈ H — the ONLY sink is Hubble-rate slow.
  (3) DRIVE: violent-relaxation deposits arrive and redistribute on the halo
      dynamical time, t_dyn = Δ_vir^{−1/2}·t_H ≈ 0.075 t_H (run_violent_relaxation).
SANDPILE MAPPING (Hwa–Kardar/BTW): non-conserved slow DRIVE is allowed (BTW's
grain addition is non-conserved too); what must be conservative is the FAST
REDISTRIBUTION, with dissipation slow/boundary-only. Here redistribution is
diffusive-conserving (1) and the only sink acts at rate H (2), so the effective
dissipation per redistribution step is
    ε = κ · t_dyn ≈ H · t_dyn = Δ_vir^{−1/2} ≈ 0.075  —  SMALL, but NOT zero.
CONSEQUENCE (mean-field SOC with dissipation ε): the avalanche law is
P(s) ∝ s^{−3/2} e^{−s/s_c} with cutoff s_c ∝ 1/ε² ≈ a few × 10² — a TRUNCATED
scale-free window of ~2.5 decades with the τ = 3/2 exponent intact. Fork B is
upgraded from "conditional on an unknown conservation law" to a QUANTITATIVE
prediction: near-critical, τ ≈ 3/2, with a computable cutoff; Fork A (no window
at all) would require ε = O(1), which the membrane coefficients do not give.

We verify the mean-field dissipative-SOC statistics with a branching process
(offspring mean 1−ε): τ → 3/2, ⟨s⟩ = 1/ε exactly, cutoff ⟨s²⟩/⟨s⟩ ∝ ε^{−2}.

HONEST FLAGS. κ ≈ H is the order-of-magnitude apparent-horizon surface gravity
(exact: κ = (1/R_AH)(1 − Ṙ_AH/2HR_AH)); the redistribution time is identified
with the deposit's own dynamical time; linearised treatment; mean-field class
justified by the all-to-all §5(ii) coupling. A full GR derivation of the
linearised height equation remains the rigorous version of step (2).
"""
import numpy as np

RNG = np.random.default_rng(5)

# ---------------------------------------------------------------------------
# 1. the horizon's dissipation parameter from membrane ingredients
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. Membrane ingredients ⇒ effective dissipation ε = κ·t_dyn = Δ_vir^(−1/2)")
print("=" * 78)
Delta_vir = 178.0
eps_horizon = Delta_vir ** -0.5
print(f"  transport: Damour–NS, dynamical-horizon ζ = +1/16πG (diffusive, locally")
print(f"             conserving) — the redistribution channel;")
print(f"  sink:      linearised Raychaudhuri relaxation at κ ≈ H (Hubble-slow);")
print(f"  drive:     violent-relaxation deposits, redistribution time t_dyn;")
print(f"  ⇒ ε = κ·t_dyn ≈ H·t_dyn = Δ_vir^(−1/2) = {eps_horizon:.3f}")
print(f"  ⇒ predicted avalanche cutoff s_c ≈ 2/ε² = {2/eps_horizon**2:.0f}")
print("  The horizon is NEARLY conservative not by assumption but because its")
print("  redistribution is dynamical-time fast while its only sink is Hubble-slow.")

# ---------------------------------------------------------------------------
# 2. mean-field dissipative SOC: branching process with offspring mean 1−ε
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. Branching-process check: τ = 3/2 window with cutoff s_c ∝ ε^(−2)")
print("=" * 78)
def avalanches(eps, n=200000, s_cap=2_000_000, gen_cap=200000):
    lam = 1.0 - eps
    active = np.ones(n, dtype=np.int64)
    s = np.ones(n, dtype=np.int64)
    for _ in range(gen_cap):
        alive = active > 0
        if not alive.any():
            break
        active[alive] = RNG.poisson(lam * active[alive])
        s += active
        np.clip(s, None, s_cap, out=s)
    return s

def tau_fit(s, s_lo, s_hi):
    """Log-binned least-squares slope of P(s) over [s_lo, s_hi]."""
    bins = np.unique(np.round(np.geomspace(s_lo, s_hi, 14)).astype(int))
    hist, edges = np.histogram(s, bins=bins)
    widths = np.diff(edges)
    centers = np.sqrt(edges[1:] * edges[:-1])
    dens = hist / widths / len(s)
    ok = dens > 0
    return -np.polyfit(np.log(centers[ok]), np.log(dens[ok]), 1)[0]

eps_list = [0.30, 0.15, eps_horizon, 0.02]
mean_s, cut_s, taus = [], [], []
print(f"    {'ε':>7} {'⟨s⟩':>9} {'1/ε':>7} {'⟨s²⟩/⟨s⟩':>10} {'2/ε²':>8} {'τ (fit)':>8}")
for eps in eps_list:
    s = avalanches(eps)
    ms = s.mean()
    cs = (s.astype(float) ** 2).mean() / ms
    t = tau_fit(s, 2, max(4, 0.5 / eps ** 2))
    mean_s.append(ms); cut_s.append(cs); taus.append(t)
    print(f"    {eps:>7.3f} {ms:>9.2f} {1/eps:>7.1f} {cs:>10.1f} {2/eps**2:>8.0f} {t:>8.2f}")
slope_cut = np.polyfit(np.log(eps_list), np.log(cut_s), 1)[0]
print(f"  ⟨s⟩·ε = {[f'{m*e:.2f}' for m, e in zip(mean_s, eps_list)]} (exact: 1)")
print(f"  cutoff scaling: ⟨s²⟩/⟨s⟩ ∝ ε^{slope_cut:.2f}   (prediction: −2)")
print(f"  τ at small ε: {taus[-1]:.2f}   (mean-field SOC: 3/2)")

# the horizon's window
s_hor = avalanches(eps_horizon)
tau_hor = tau_fit(s_hor, 2, 0.5 / eps_horizon ** 2)
window = np.log10(2 / eps_horizon ** 2)
print(f"\n  AT THE HORIZON's ε = {eps_horizon:.3f}:  τ = {tau_hor:.2f} over a "
      f"~{window:.1f}-decade window, cutoff s_c ≈ {2/eps_horizon**2:.0f}")

# ---------------------------------------------------------------------------
# 3. the Fork A/B dichotomy, now quantitative
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. Fork A vs Fork B becomes a NUMBER, not an unknown")
print("=" * 78)
s_a = avalanches(0.5)         # Fork-A-like: order-unity dissipation
tau_a = tau_fit(s_a, 2, 30)
print(f"  ε = O(1) (Fork A, no conservation): ⟨s⟩ = {s_a.mean():.1f}, no scale-free")
print(f"    window (effective slope {tau_a:.2f} over a decade — dominated by cutoff);")
print(f"  ε = {eps_horizon:.3f} (membrane value):  τ ≈ 3/2 over ~{window:.1f} decades, then cutoff;")
print(f"  ε → 0 (pure SOC): unbounded window (verified: cutoff ∝ ε^{slope_cut:.1f}).")
print("  ⇒ the membrane coefficients place the horizon in the TRUNCATED-critical")
print("    regime: the τ ≈ 3/2 signature at z* survives, but with a predicted")
print("    cutoff — a SHARPER, more falsifiable Fork B than the unconditional one.")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — the OPEN conservation-law row can be closed (prototype level)")
print("=" * 78)
print(f"""  The ledger's 'conservation law: OPEN (decides Fork A vs B)' is answerable
  from horizon hydrodynamics the paper already cites (§7.2 membrane paradigm):
    • redistribution is locally conserving (Damour–NS, ζ_AH = +1/16πG > 0);
    • the only sink is Hubble-rate slow (linearised Raychaudhuri, κ ≈ H);
    • hence ε = κ·t_dyn ≈ Δ_vir^(−1/2) = {eps_horizon:.3f} — NEARLY conservative,
      automatically, with a TRUNCATED scale-free window: τ ≈ 3/2 (measured
      {tau_hor:.2f}) up to s_c ≈ 2/ε² ≈ {2/eps_horizon**2:.0f}, i.e. ~{window:.1f} decades.
  Fork B is thus not 'conditional on an unknown law' but PREDICTED in truncated
  form; pure Fork A (no window) would need ε = O(1), which the membrane
  coefficients do not give. Manuscript action: replace the OPEN row with
  'NEARLY-CONSERVATIVE (ε = H·t_dyn ≈ 0.075, membrane-derived) ⇒ truncated
  τ = 3/2 window, cutoff s_c ~ 2/ε²' — and flag the full linearised-GR
  derivation of the height equation as the remaining rigorous step.""")

# validation
for m, e in zip(mean_s, eps_list):
    assert abs(m * e - 1.0) < 0.1, f"mean avalanche size must be 1/eps (eps={e})"
assert -2.5 < slope_cut < -1.5, "cutoff must scale as eps^-2 (mean-field SOC)"
assert 1.3 < taus[-1] < 1.7, "small-eps exponent must be the mean-field tau=3/2"
assert 1.3 < tau_hor < 1.8, "horizon-eps window must still show tau~3/2"
assert eps_horizon < 0.1, "membrane dissipation parameter must be small (near-conservative)"
print("\n[validate] dns-roughening assertions passed.")
