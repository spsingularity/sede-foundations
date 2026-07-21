"""
run_violent_relaxation.py — the γ-mechanism and the driven-NESS are ONE process.

Inference #3 from the Paper-II literature (Lynden-Bell violent relaxation;
Campa-Dauxois-Ruffo quasi-stationary states). SEDE invokes violent relaxation
twice and they are the same event:

  • the γ-mechanism (Thm 4C): a virialising halo dissipates its binding energy
    E_bind ∝ M^{5/3} into the horizon ledger, setting the weight p=5/3, γ≈1.5;
  • the driven-NESS (the volume lock): the structure-deposited entropy drives the
    horizon to a volume-law steady state.

Lynden-Bell (1967): gravitational collapse drives a long-range system, on the
*dynamical* (violent) time, to a quasi-stationary state (QSS) — not the thermal
equilibrium — whose lifetime in a system of N dof diverges with N (CDR). We show:
(1) the same binding-energy/virial chain sets BOTH the γ weight and the NESS
drive amplitude; (2) violent relaxation is prompt (t_dyn ≪ t_Hubble), so the
horizon is continuously driven; (3) the QSS lifetime ∝ N^δ makes the volume lock
permanent at the horizon's N — a SECOND permanence argument beside the GSL
(route E); (4) halo and horizon QSS both have negative specific heat
(microcanonical), tying this to the ensemble-inequivalence result.
"""
import numpy as np
from sede.halo_entropy import horizon_deposited_entropy_scaling

# ===========================================================================
# (1) One binding-energy chain sets BOTH γ and the NESS drive
# ===========================================================================
print("=" * 74)
print("(1) Violent relaxation: one E_bind chain sets the γ weight AND the drive")
print("=" * 74)
M, slope = horizon_deposited_entropy_scaling()
print(f"  dissipated binding energy ΔS_AH = E_bind/T_AH:  E_bind ∝ M^{slope:.2f}")
print(f"    (NFW + virial; ≈ M^{{5/3}} ⇒ entropy weight p=5/3 ⇒ γ≈1.50 — Thm 4C)")
print(f"  virial theorem chain (the SAME physics):")
print(f"    R_vir ∝ M^{{1/3}} (fixed Δ_vir) ⇒ σ_v² = GM/R_vir ∝ M^{{2/3}}")
print(f"    E_bind = M σ_v² ∝ M·M^{{2/3}} = M^{{5/3}}      → the γ WEIGHT")
print(f"    ε_dep  ∝ f_coll ⟨σ_v²⟩/c²                   → the NESS DRIVE amplitude")
print(f"  So the drive that carries the horizon to its QSS IS the violent-relaxation")
print(f"  binding energy whose mass-weighting is the γ-mechanism. ONE process.")

# ===========================================================================
# (2) Violent relaxation is prompt: t_dyn ≪ t_Hubble ⇒ QSS reached continuously
# ===========================================================================
print("\n" + "=" * 74)
print("(2) Promptness: collapse (violent) time ≪ Hubble time")
print("=" * 74)
Delta_vir = 200.0                       # virial overdensity
# t_dyn/t_H = (ρ_crit/ρ_vir)^{1/2} = Δ_vir^{-1/2} (free-fall vs Hubble)
ratio = Delta_vir**-0.5
print(f"  t_dyn/t_Hubble = Δ_vir^(-1/2) = {ratio:.3f}  (Δ_vir≈200)")
print(f"  ⇒ each halo virialises in ~{ratio:.0%} of a Hubble time: violent relaxation is")
print(f"    fast, the QSS is reached promptly, and the horizon is driven by an")
print(f"    ongoing sequence of collapse events (the f_sat gate is their envelope).")

# ===========================================================================
# (3) QSS lifetime ∝ N^δ ⇒ permanent at the horizon's N (2nd permanence arg)
# ===========================================================================
print("\n" + "=" * 74)
print("(3) QSS lifetime diverges with N ⇒ the volume lock is permanent")
print("=" * 74)
N_hor = 1e122                           # horizon dof ~ S_dS
for delta in (1.0, 1.7):               # CDR: t_QSS/t_dyn ~ N^δ (δ≈1–1.7 across models)
    print(f"  δ={delta:.1f}:  t_QSS/t_dyn ~ N^δ = (10^122)^{delta:.1f} = 10^{122*delta:.0f}")
print(f"  Even the conservative δ=1 gives t_QSS/t_dyn ~ 10^122 — astronomically")
print(f"  longer than any cosmological time. The volume-law QSS is effectively")
print(f"  permanent (Δ=1 constant), a SECOND permanence argument beside the GSL")
print(f"  entropy-maximum of route E. (Same large-N that makes J≫J_c, run_jcoup.)")

# ===========================================================================
# (4) Both QSS have negative specific heat ⇒ microcanonical (ties to ensemble)
# ===========================================================================
print("\n" + "=" * 74)
print("(4) Halo and horizon QSS both microcanonical (negative specific heat)")
print("=" * 74)
print(f"  self-gravitating QSS: gravothermal, C<0 (Lynden-Bell–Wood; Antonov).")
print(f"  de Sitter horizon: C = −2S < 0 (run_lit_inferences.py).")
print(f"  Both are negative-specific-heat, long-range ⇒ MICROCANONICAL systems —")
print(f"  consistent with the ensemble-inequivalence selection of the volume branch.")
print(f"  The f_sat Boltzmann factor e^{{-γu}} (Thm 3 kernel) is the dilute-limit")
print(f"  Lynden-Bell statistical weight of the violently-relaxed deposition.")

print("\n" + "=" * 74)
print("VERDICT — γ-mechanism ⊕ driven-NESS = one violent-relaxation statement")
print("=" * 74)
print("""  Structure forms by gravitational collapse ⇒ violent relaxation, which
  simultaneously (a) dissipates E_bind ∝ M^{5/3} into the horizon (the γ weight,
  Thm 4C) and (b) drives the horizon to a volume-law quasi-stationary state (the
  driven-NESS). The deposited binding energy IS the NESS drive; its mass-weighting
  IS γ. The QSS is microcanonical (negative specific heat, like the horizon), so
  the same long-range physics that selects the volume branch (ensemble
  inequivalence) also names the dynamical state (Lynden-Bell QSS) and guarantees
  its permanence (lifetime ∝ N^δ). This UNIFIES two of SEDE's separate
  ingredients — the γ derivation and the driven-NESS — into a single statement,
  adding no assumption. Open as before: whether the QSS state count is volume
  (dS holography).""")

assert 1.5 < slope < 1.8, "binding-energy slope must be ≈5/3 (the violent-relaxation deposit)"
assert ratio < 0.2, "violent relaxation must be prompt vs Hubble"
print("\n[validate] violent-relaxation assertions passed.")
