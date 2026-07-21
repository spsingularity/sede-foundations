"""
run_algebra_type.py — why "everything equilibrium points to area" is an algebra-
type statement, not a law of horizons: the area CEILING is the Type II₁
max-entropy property of the ETERNAL static patch, which a driven FRW horizon
need not share.

THE GAP. §7/§7.1 concede that the equilibrium handles (DSSYK in the Narovlansky–
Verlinde dictionary, the dS modular Hamiltonian, black-hole thermodynamics) all
favour the AREA count, and defuse this with a state-dependent-count argument the
paper itself calls "an effective-model statement one layer out."

THE CLOSURE (algebraic sharpening, CLPW = Chandrasekaran–Longo–Penington–Witten
2023). In the crossed-product construction:
  • the ETERNAL dS static patch (+ observer) has a Type II₁ algebra — the unique
    case possessing a MAXIMUM-entropy state, and that maximum IS the
    Gibbons–Hawking area. The area bound is the II₁ ceiling.
  • the black-hole exterior gives Type II_∞ — NO maximum-entropy state; entropy
    unbounded above.
So the "S ≤ Area" ceiling is a special property of the equilibrium (eternal,
matter-free) algebra, not a universal horizon law. A matter-filled FRW apparent
horizon is NOT the eternal static patch; whether its algebra is II₁ or II_∞ is a
sharply-posed open question — and if II_∞, there is no ceiling for the driven
horizon to violate. DSSYK inherits the same status: its Hilbert space is 2^{N/2}
BY CONSTRUCTION (S_max = (N/2)ln2), so it cannot pose super-extensivity; whether
N itself scales as R² or R³ is a UV input to the dictionary, not an SYK output.

TOY DEMONSTRATION (finite-dimensional analog — a tracial algebra on a FIXED
Hilbert space has max entropy ln(dim), reached by the tracial state; ENLARGING
the algebra (a horizon absorbing structure dof) removes any fixed ceiling):
  (1) fixed register: Page-random states saturate near, never exceed, ln(dim);
  (2) growing register: the ceiling itself grows linearly — no maximum state.
This is the II₁-vs-II_∞ dichotomy in its finite-dimensional shadow.

HONEST FLAG. The toy illustrates the dichotomy; it does not compute the FRW
horizon's algebra type. The manuscript gain is the REFRAMING: "equilibrium ⇒
area" becomes "the equilibrium algebra is the one with a ceiling" — a precise,
citable statement (CLPW) replacing the effective-model hedge, plus one new
sharply-posed question (the FRW patch's type) for the residue ledger.
"""
import numpy as np

RNG = np.random.default_rng(2)

def page_entropy(nA, nB, trials=30):
    """Mean entanglement entropy of Haar-random pure states, subsystem A."""
    dA, dB = 2 ** nA, 2 ** nB
    S = []
    for _ in range(trials):
        psi = RNG.standard_normal((dA, dB)) + 1j * RNG.standard_normal((dA, dB))
        psi /= np.linalg.norm(psi)
        lam = np.linalg.svd(psi, compute_uv=False) ** 2
        lam = lam[lam > 1e-15]
        S.append(float(-(lam * np.log(lam)).sum()))
    return float(np.mean(S))

# ---------------------------------------------------------------------------
# 1. fixed algebra (II₁ analog): a ceiling exists and Page states saturate it
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. FIXED register (II₁ analog): max entropy = ln(dim), Page states saturate")
print("=" * 78)
n = 8
ceiling = n * np.log(2)
S_page = page_entropy(n, n)
S_page_theory = n * np.log(2) - 0.5          # Page: ln dA − dA/2dB, dA=dB
print(f"  horizon register n = {n} qubits: ceiling S_max = n ln2 = {ceiling:.3f}")
print(f"  Haar-random states: ⟨S⟩ = {S_page:.3f}   (Page value ln dA − 1/2 = {S_page_theory:.3f})")
print(f"  maximally mixed (tracial) state: S = ln(dim) = {ceiling:.3f} — the CEILING")
print("  ⇒ a fixed (equilibrium) algebra has a maximum-entropy state; for the")
print("    eternal dS static patch that ceiling IS the Gibbons–Hawking area (II₁).")

# ---------------------------------------------------------------------------
# 2. growing algebra (II_∞ analog): the ceiling is unbounded
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("2. GROWING register (II_∞ analog): absorbing structure dof removes the ceiling")
print("=" * 78)
S_seq, ceilings = [], []
print(f"    {'step':>5} {'register n(t)':>14} {'ceiling n ln2':>14} {'⟨S⟩ random':>11}")
for t, nt in enumerate([4, 6, 8, 10]):
    c = nt * np.log(2)
    Sp = page_entropy(nt, nt, trials=15)
    S_seq.append(Sp); ceilings.append(c)
    print(f"    {t:>5} {nt:>14} {c:>14.3f} {Sp:>11.3f}")
grow = np.polyfit(range(4), ceilings, 1)[0]
print(f"  ceiling growth per absorption step: {grow:.3f} nats — UNBOUNDED; no")
print("  maximum-entropy state exists for the growing (driven) system.")
print("  ⇒ a horizon that absorbs structure dof is not described by the fixed")
print("    equilibrium algebra; the II₁ area ceiling simply does not apply to it.")

# ---------------------------------------------------------------------------
# 3. the DSSYK corollary
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("3. DSSYK corollary: 'DSSYK ⇒ area' is built into the frame")
print("=" * 78)
for N_syk in (100, 400, 1600):
    print(f"    N = {N_syk:>5} Majoranas: dim H = 2^(N/2), S_max = (N/2)ln2 = "
          f"{0.5 * N_syk * np.log(2):>7.1f} — extensive in N by construction")
print("  ⇒ DSSYK cannot exhibit super-extensive (volume, N^{3/2}) entropy for ANY")
print("    state: its Hilbert space is fixed at 2^{N/2}. Whether the boundary N")
print("    scales as R² (area) or R³ (volume) is a UV input to the dictionary,")
print("    not an SYK output — so 'DSSYK favours area' has no force against a")
print("    driven horizon; it restates the frame's own counting.")

# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("VERDICT — the equilibrium-area tension, restated as an algebra type")
print("=" * 78)
print("""  The standing objection 'every concrete handle gives area' is the statement
  'the ETERNAL static patch's algebra is Type II₁, whose maximum entropy is the
  Gibbons–Hawking area' (CLPW 2023). That ceiling is a property of the
  equilibrium algebra, not of horizons per se: the black-hole exterior is
  already II_∞ (no ceiling), and the matter-filled FRW apparent horizon is not
  the eternal patch. The manuscript's state-dependent count (§7.1) can therefore
  be restated precisely: equilibrium horizon ⇔ fixed II₁ algebra ⇔ area ceiling;
  driven horizon ⇔ enlarging algebra ⇔ no ceiling — with ONE new sharply-posed
  question for the ledger: the TYPE of the FRW static-patch algebra (II₁ or
  II_∞). This replaces the 'effective-model, one layer out' hedge with a citable
  algebraic dichotomy and localises what remains open.
  FLAG: the finite-dimensional toy ILLUSTRATES the dichotomy; computing the FRW
  algebra type is the (well-posed) open problem handed to dS holography.""")

# validation
assert S_page < ceiling, "random-state entropy must stay below the ceiling"
assert abs(S_page - S_page_theory) < 0.1, "Page saturation on the fixed register"
assert grow > 1.0, "growing register's ceiling must grow (no maximum state)"
assert all(S_seq[i] < ceilings[i] for i in range(4)), "S below the running ceiling"
print("\n[validate] algebra-type assertions passed.")
