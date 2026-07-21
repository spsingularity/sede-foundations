"""
run_pbh_driving.py — does ANY black hole get driven enough to lock volume?
The serious test of route 1 (the would-be "Paper III" hook).

The most-driven BH candidate is not an accreting BH but a PRIMORDIAL BH AT
FORMATION: a PBH forms when a horizon-scale overdensity collapses, so — unlike an
existing BH that accretes a small fraction — the ENTIRE horizon is built by a
violent-relaxation event. If any BH activates volume, this is where.

The drive on the area↔volume order parameter is the structure entropy DEPOSITED
relative to the horizon entropy created, and it must reach the bistability barrier
h_c ~ O(1) (even after near-critical amplification) to lock volume. At PBH
formation the deposited entropy is the radiation entropy of the collapsing Hubble
patch; the created entropy is the BH entropy:

    S_BH  ~ (M/M_P)^2          (Bekenstein–Hawking)
    S_rad ~ (M/M_P)^{3/2}      (radiation entropy of the horizon patch at formation)
    f_form = S_rad / S_BH ~ (M/M_P)^{-1/2}

so the formation drive falls as (M_P/M)^{1/2} — tiny for any super-Planckian PBH.
We check it against the barrier, with and without the cosmic near-critical
amplification, across the PBH mass range.
"""
import numpy as np
from sede import barrow_bh as bh

M_P_g = 2.176e-5                                   # Planck mass (g)
h_c = 0.42                                         # bistability barrier (run_deposition_drive)
chi = 6.66e6                                       # near-critical amplification (CHR identity)

print("=" * 74)
print("Route 1 test — PBH formation drive vs the volume-locking barrier")
print("=" * 74)
print(f"  barrier h_c = {h_c};  near-critical amplification χ ≈ {chi:.1e}")
print(f"  {'PBH mass':>12}{'M/M_P':>11}{'f_form=(M_P/M)^½':>18}{'amplified χ·f':>14}  locks?")
masses_g = [5e14, 1e20, 1e30, 1.989e33, 1.989e41]   # ~evap PBH, asteroid, Earth-ish, Msun, SMBH
labels   = ["evap PBH ~5e14 g", "1e20 g", "1e30 g", "1 Msun", "1e8 Msun SMBH"]
for Mg, lab in zip(masses_g, labels):
    x = Mg / M_P_g
    f_form = x ** -0.5                              # S_rad/S_BH ~ (M/M_P)^{-1/2}
    amp = chi * f_form
    locks = (amp > h_c)
    print(f"  {lab:>12}{x:>11.1e}{f_form:>18.1e}{amp:>14.1e}  {'YES' if locks else 'no (area)'}")
print("  even at PBH formation — the maximally-driven case — the deposited structure")
print("  entropy is ≪ the BH entropy created (S_rad/S_BH ~ (M/M_P)^{-1/2}), so the")
print("  drive falls vastly short of the barrier for every super-Planckian BH, even")
print("  with the full near-critical amplification and even as a one-shot event.")

# the cosmic horizon, for contrast, integrates a SUSTAINED drive to f~1
print("\n  contrast — cosmic apparent horizon: a SUSTAINED structure drive integrated")
print("  over a Hubble time reaches cumulative f_sat → 1 (run_deposition_drive.py),")
print("  clearing h_c at z*. A PBH gets ONE tiny formation kick. No comparison.")

print("\n" + "=" * 74)
print("VERDICT — is there a non-null black-hole prediction? (the Paper-III question)")
print("=" * 74)
print("""  NO. The most-driven black hole conceivable — a PBH built entirely by a
  collapsing overdensity — still forms AREA-law, because the structure entropy
  deposited at formation is negligible against the black-hole entropy created
  (f_form ~ (M_P/M)^{1/2} ≪ h_c even amplified). Accreting/merging BHs are weaker
  still. So:

  • SEDE robustly predicts ALL black holes are area-law (Δ_BH = 0) ⇒ standard
    Hawking evaporation — the earlier prediction holds, now from a calculation,
    not a hand-wave.
  • The cosmic apparent horizon is the UNIQUE volume-law horizon, because it is
    the only one driven *sustainedly* (the whole structure-formation history) to
    f_sat → 1. Its uniqueness is a result, not an assumption.
  • Therefore there is NO non-null, BH-specific observable — the would-be Paper III
    has no hook. The BH content is exactly: the null prediction (BHs are area) plus
    the falsifier (a volume-law isolated BH would refute the state-dependent
    mechanism). That belongs in Paper II's §7.1, not a separate paper.""")

assert (chi * (5e14 / M_P_g) ** -0.5) < h_c, "even amplified, the evap-PBH formation drive < barrier"
assert (5e14 / M_P_g) ** -0.5 < 1e-8, "formation drive must be tiny for a 5e14 g PBH"
print("\n[validate] PBH-driving assertions passed.")
