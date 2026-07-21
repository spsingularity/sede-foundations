"""
run_chr_amplitude.py — the deferred CHR amplitude, computed: the z*-localised
signal has TWO layers, and the membrane number ε = H·t_dyn ≈ 0.075 sets the
detectable one.

THE GAP. §6 defers "a quantitative amplitude and survey signal-to-noise" of the
z*-localised excess to the companion CHR analysis — the last un-forecast piece of
the prediction stack.

THE CLOSURE (both layers, honestly separated):
  (1) DETERMINISTIC LAG (detectable). Paper I's f_sat is the ADIABATIC limit of a
      relaxational order-parameter dynamics (sede/chr_mechanism.py). The
      relaxation rate is no longer free: deposits redistribute on the dynamical
      time, so Γ = 1/(H·t_dyn) = 1/ε ≈ 13.3 per e-fold — the SAME membrane number
      behind the SOC cutoff (run_dns_roughening.py). The first-order effect is a
      LAG of the gate behind its equilibrium, Δφ ≈ ε·dφ_eq/dn, peaking where the
      gate moves fastest (the z* epoch) with amplitude ∼ ε·max(dφ_eq/dn) — a few
      per cent. Propagated through the λ = ½ closure to H(z) and the growth rate,
      this is a percent-level, z*-localised, growth-locked distortion of fσ₈(z):
      forecastable against Euclid/DESI-class precision (σ ≈ 1% per z-bin).
  (2) STOCHASTIC SCALE-FREE LAYER (NOT detectable at 2-point level — honest
      negative). The avalanche shot noise carries the τ ≈ 3/2 statistics but rms
      fractional amplitude ∼ ε_dep·√s_c ∼ 10⁻⁶ — six orders below survey
      precision. The τ = 3/2 exponent is therefore NOT accessible through mean
      background/growth measurements; it would require the fluctuation layer
      (higher-order/ISW statistics) — we pre-register that limitation instead of
      letting P3–P5 imply otherwise.
NET: the §6 amplitude deferral is closed: the detectable CHR observable is the
deterministic ε-lag (percent level, S/N computed below); the scale-free exponent
rides on a 10⁻⁶ layer and is out of reach of two-point survey statistics.

HONEST FLAGS. Whether a deposit "counts" in f_sat at deposition or after
redistribution is an O(ε) modelling choice — the lag amplitude is the size of
that choice, so this forecast is the honest uncertainty band of the gate as much
as a new signal; the background closure used here is the algebraic λ = ½ fixed
point (not the full ODE background); survey σ is a round 1% per bin.
"""
import numpy as np
from sede.chr_mechanism import (f_eq, deposited_entropy_fraction, transition_redshift,
                                GAMMA)
from sede.friedmann import compute_growth_factor, compute_growth_model

OM, OR = 0.30, 9e-5
ODE = 1.0 - OM - OR
EPS = 178.0 ** -0.5                     # membrane dissipation ε = H·t_dyn
GAM_RELAX = 1.0 / EPS                   # relaxation rate per e-fold (fixed, not free)

# ---------------------------------------------------------------------------
# 1. the lagged gate: dφ/dn = Γ(φ_eq − φ) with Γ = 1/ε (membrane-fixed)
# ---------------------------------------------------------------------------
print("=" * 78)
print(f"1. The lagged gate: Γ = 1/ε = {GAM_RELAX:.1f}/e-fold (membrane-fixed, not tuned)")
print("=" * 78)
n = np.linspace(-4.0, 0.0, 4001)                     # n = ln a
z_of_n = np.exp(-n) - 1.0
D = compute_growth_factor(z_of_n, OM)
x = D ** 2
phi_eq = f_eq(x)
phi = np.empty_like(n); phi[0] = phi_eq[0]
dn = n[1] - n[0]
for k in range(1, len(n)):
    phi[k] = phi[k - 1] + dn * GAM_RELAX * (phi_eq[k - 1] - phi[k - 1])
dphi = phi_eq - phi                                   # the lag (>0 while rising)
dphidn = np.gradient(phi_eq, n)
k_pk = int(np.argmax(dphi))
z_star = transition_redshift(OM)
print(f"  lag Δφ = φ_eq − φ_lag peaks at z = {z_of_n[k_pk]:.2f}  (z* = {z_star:.2f})")
print(f"  peak amplitude Δφ = {dphi[k_pk]:.4f}   first-order estimate ε·max(dφ/dn) = "
      f"{EPS * dphidn.max():.4f}")
print(f"  fractional: Δφ/φ at peak = {dphi[k_pk]/max(phi_eq[k_pk],1e-9):.3f}")
print("  ⇒ a few-percent, z*-localised distortion of the gate — the size of the")
print("    'deposit counts now vs after redistribution' choice, i.e. the honest")
print("    O(ε) uncertainty AND the CHR signal, one number.")

# ---------------------------------------------------------------------------
# 2. propagate to H(z) and fσ8(z); Euclid/DESI-class S/N
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. Background + growth propagation, and the survey forecast")
print("=" * 78)
def make_E(phi_grid):
    """λ=1/2 algebraic closure: E² = Ω_m(1+z)³ + Ω_r(1+z)⁴ + c·φ·E, E(0)=1."""
    c = (1.0 - OM - OR) / max(phi_grid[-1], 1e-9)     # calibrate at z=0
    zg = z_of_n[::-1]; pg = phi_grid[::-1]            # ascending z for interp
    def E_of(z):
        z = np.atleast_1d(np.asarray(z, float))
        p = np.interp(z, zg, pg, right=0.0)
        m = OM * (1 + z) ** 3 + OR * (1 + z) ** 4
        return 0.5 * (c * p + np.sqrt((c * p) ** 2 + 4 * m))
    return E_of

E_ad, E_lag = make_E(phi_eq), make_E(phi)
z_q = np.linspace(0.0, 3.0, 61)
dH = E_lag(z_q) / E_ad(z_q) - 1.0
k_h = int(np.argmax(np.abs(dH)))
print(f"  |δH/H| peaks at z = {z_q[k_h]:.2f}: {abs(dH[k_h])*100:.2f}%   (BAO ~0.5%/bin)")

D_a, g_a = compute_growth_model(z_q, OM, E_ad)
D_l, g_l = compute_growth_model(z_q, OM, E_lag)
dfs8 = (g_l * D_l) / (g_a * D_a) - 1.0                # fractional δ(fσ8)
zb = np.arange(0.3, 1.75, 0.2)                        # Euclid-like bins
db = np.interp(zb, z_q, dfs8)
SIG = 0.01                                            # 1% per bin
snr = float(np.sqrt(np.sum((db / SIG) ** 2)))
print(f"    {'z bin':>6} {'δfσ8/fσ8':>10} {'S/N':>6}")
for zz, dd in zip(zb, db):
    print(f"    {zz:>6.1f} {dd:>10.4f} {abs(dd)/SIG:>6.1f}")
print(f"  total S/N (8 bins, σ = 1%): {snr:.1f}")
print("  ⇒ the deterministic ε-lag is a percent-level, z*-localised, growth-locked")
print("    signature — at the edge of Euclid/DESI-class sensitivity, with its")
print("    amplitude set by the SAME ε as the SOC cutoff (one number, two roles).")

# ---------------------------------------------------------------------------
# 3. the stochastic layer: honest negative for the τ = 3/2 exponent at 2-point
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. The scale-free layer: amplitude ε_dep·√s_c — out of reach (pre-registered)")
print("=" * 78)
eps_dep = deposited_entropy_fraction(OM)
s_c = 2.0 / EPS ** 2
rms = eps_dep * np.sqrt(s_c)
print(f"  ε_dep = {eps_dep:.2e},  s_c = 2/ε² = {s_c:.0f}  ⇒  rms fractional shot noise")
print(f"  ~ ε_dep·√s_c = {rms:.1e} — six orders below the 1% survey precision.")
print("  ⇒ HONEST NEGATIVE: the τ ≈ 3/2 exponent does NOT imprint on mean")
print("    background/growth two-point statistics at any plausible precision; it")
print("    lives in the fluctuation layer only. P3–P5 should be phrased as the")
print("    ε-lag shape (detectable) + the scale-free statistics (theoretical")
print("    discriminator, not a near-term survey observable).")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — the §6 amplitude deferral is closed, both ways")
print("=" * 78)
print(f"""  Detectable: the deterministic lag of the gate behind its equilibrium,
  amplitude ε·max(dφ_eq/dn) ≈ {EPS*dphidn.max():.3f} at z ≈ {z_of_n[k_pk]:.1f} (≈ z*), propagating to a
  {abs(dH[k_h])*100:.1f}% peak in H(z) and Δ(fσ8) with total S/N ≈ {snr:.0f} over Euclid-like bins —
  with NO new parameter: Γ = 1/ε is the membrane number of run_dns_roughening.
  Not detectable (pre-registered): the τ ≈ 3/2 scale-free statistics, riding on a
  10⁻⁶ shot-noise layer. Manuscript action: replace §6's 'amplitude ... deferred'
  with these two numbers, and re-scope P3–P5 accordingly.""")

# validation
assert abs(eps_dep - 1.5e-7) < 5e-8, "CHR deposited-entropy fraction ~1.5e-7"
assert 0.005 < dphi[k_pk] < 0.1, "lag amplitude must be percent-level"
assert abs(z_of_n[k_pk] - z_star) < 0.8, "lag must peak near z*"
assert dphi[k_pk] < 2.5 * EPS * dphidn.max() and dphi[k_pk] > 0.4 * EPS * dphidn.max(), \
    "lag must match the first-order eps*dphi/dn estimate"
assert snr > 1.0, "the deterministic layer must be at least marginally detectable"
assert rms < 1e-4, "the stochastic layer must be far below survey precision"
print("\n[validate] chr-amplitude assertions passed.")
