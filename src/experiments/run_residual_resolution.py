"""
run_residual_resolution.py — how to resolve the round-2 residual.

THE RESIDUAL (objection, round 2): the area-law-theorem defence of V1 is a
*conditional* — it treats the horizon as a many-body system whose dof + Hamiltonian
are not actually defined; defining them is the open dS-holography problem, and the
standard equilibrium handles (the de Sitter modular Hamiltonian = a LOCAL boost;
DSSYK in the Narovlansky–Verlinde dictionary; CKN; black-hole thermodynamics) all
point to AREA. So "everything points to area" looks like a standing tension.

THE RESOLUTION — the count is STATE-DEPENDENT, not a universal constant of QG.
The very discriminator the paper already uses (equilibrium → area, driven → volume)
says the area results are not counter-evidence: they are all EQUILIBRIUM horizons,
correctly area. The cosmic horizon is the one DRIVEN horizon, and the f_sat gate is
precisely the fraction of the volume capacity that structure-driving activates:

      effective count:  N_eff = N_area + (N_vol − N_area)·f_sat
      equilibrium (no structure, f_sat→0)  ⇒  area   (BH, DSSYK, modular Ham, CKN)
      driven  (structure, f_sat→1, LOCKED) ⇒  volume (the cosmic horizon)

This is NOT the rejected running-Δ background: the cosmic horizon is driven past the
spinodal early and LOCKS at Δ=1 (constant-Δ background; socatt/violrel permanence).
The state-dependence is ACROSS horizons (driven cosmic vs equilibrium BH), not a
running Δ in the cosmic background.

CONSEQUENCE — the residual becomes DOUBLY falsifiable, on two independent horizons:
  (1) cosmic Δ = 1  (driven)      — DESI DR3 + Euclid (the EOS/growth);
  (2) BH    Δ = 0  (equilibrium)  — PBH evaporation: SEDE PREDICTS PBHs evaporate
      NORMALLY (standard Hawking), because black holes are undriven.
A universal-volume-QG theory predicts the opposite for (2) — PBHs never evaporate.
So the PBH channel is not a threat to SEDE; it is a test of its MECHANISM, with a
definite prediction that DISTINGUISHES driven-NESS-SEDE from universal-volume-QG.
"""
import numpy as np
from sede import barrow_bh as bh

Msun_g = 1.989e33

# ---------------------------------------------------------------------------
# (1) The state-dependent count dissolves the "everything points to area" tension
# ---------------------------------------------------------------------------
print("=" * 74)
print("(1) The count is STATE-DEPENDENT: equilibrium→area, driven→volume")
print("=" * 74)
print("  N_eff = N_area + (N_vol − N_area)·f_sat   (f_sat = activated volume fraction)")
for label, f in [("black hole / DSSYK / dS modular Ham (equilibrium)", 0.0),
                 ("cosmic horizon today (structure-driven, locked)", 1.0)]:
    Delta_eff = f                                   # area Δ=0 (f=0) → volume Δ=1 (f=1)
    print(f"    {label:<48} f_sat={f:.0f} ⇒ Δ_eff={Delta_eff:.0f} "
          f"({'AREA' if f < 0.5 else 'VOLUME'})")
print("  ⇒ the area results (DSSYK, modular Hamiltonian, CKN, BH) are all EQUILIBRIUM")
print("    horizons — correctly area, NOT counter-evidence. Only the driven cosmic")
print("    horizon activates the volume count. The 'tension' dissolves.")

# ---------------------------------------------------------------------------
# (2) The PBH discriminator: SEDE predicts BH Δ=0 ⇒ normal evaporation
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("(2) PBH evaporation — a test of the MECHANISM (driven-NESS vs universal QG)")
print("=" * 74)
M_star_g = 5.1e14                                   # PBH mass evaporating ~now (Hawking)
M = M_star_g / Msun_g
Tr1 = bh.bh_temperature_ratio(M, 1.0)
tr1 = bh.evaporation_time_ratio_log10(M, 1.0)
print(f"  reference PBH mass M* ≈ {M_star_g:.1e} g (standard Hawking: evaporates within ~1/H₀)")
print(f"  {'hypothesis':<34}{'Δ_BH':>5}{'T_B/T_H':>12}{'t_evap ratio':>14}   PBH evaporates?")
print(f"  {'driven-NESS SEDE (BH undriven)':<34}{0:>5}{1.0:>12.1e}{0.0:>14.0f}   YES (standard)")
print(f"  {'universal-volume QG':<34}{1:>5}{Tr1:>12.1e}{('×10^%.0f'%tr1):>14}   NO (never)")
print("  ⇒ SEDE-driven-NESS PREDICTS standard PBH evaporation (Δ_BH=0); a detection of")
print("    the ~MeV Hawking γ-ray signature CONFIRMS the equilibrium-area side and")
print("    DISTINGUISHES SEDE from a universal-volume theory (which forbids it).")
print("    (The black hole being area-law is thus a SEDE prediction, not a problem.)")

# ---------------------------------------------------------------------------
# (3) The two-channel empirical resolution + the residual theoretical step
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("VERDICT — how the residual is resolved")
print("=" * 74)
print("""  EMPIRICALLY (available now, no dS holography needed): the count is pinned on
  TWO independent horizons —
    • cosmic Δ = 1  (driven)     ← DESI DR3 + Euclid, σ(Δ)~0.09 (Paper I §6);
    • BH    Δ = 0  (equilibrium) ← PBH evaporation (standard Hawking γ-rays).
  Confirming both establishes the *state-dependent* count as fact AND confirms the
  equilibrium-vs-driven mechanism — and removes the 'DSSYK/modular-Ham/BH point to
  area' tension, since those are exactly the equilibrium (Δ=0) horizons SEDE expects.

  THEORETICALLY (the honest remaining step): prove that a structure-DRIVEN horizon
  has a volume-law steady-state count while the equilibrium one is area — i.e. the
  driven-NESS conjecture, now sharpened to 'the f_sat gate is the activated-volume
  fraction.' This still needs a specified horizon-dof model (the conditional of V1),
  but it is no longer in tension with the equilibrium area results; they are its
  f_sat→0 limit. The residual is thus reduced from 'derive the count' to 'prove the
  driven horizon activates the volume capacity' — with the empirical two-channel
  test deciding it either way regardless.""")

# validation
assert bh.evaporation_time_ratio_log10(M, 1.0) > 30, "Δ=1 must forbid PBH evaporation (huge t ratio)"
assert abs(bh.evaporation_time_ratio_log10(M, 0.0)) < 1e-9, "Δ=0 = standard Hawking (ratio 1)"
print("\n[validate] residual-resolution assertions passed.")
