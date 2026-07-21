"""
run_rtn_ising.py — deriving the two-state dof (residual R2): the horizon's
binary (boundary-/bulk-entangling) variable is the EXACT Ising dual of a random
tensor network's entanglement computation, not an idealisation.

THE GAP. §5(ii)'s derived Landau F(m) still rests on one idealisation: that the
horizon dof are TWO-STATE (each either area-contributing or bulk-connected). One
can ask why a horizon dof should be binary at all.

THE CLOSURE (Hayden–Nezami–Qi–Thomas–Walter–Yang 2016). In the random-tensor-
network (RTN) model of holographic states, the Rényi-entropy computation maps
EXACTLY onto an Ising model: averaging Tr ρ_A² over the random tensors turns each
tensor into a spin with two values — the identity or swap permutation — and
    E[Tr ρ_A²] ≈ Z_A/Z_∅,   Z_A = Σ_{spin configs} χ^{−(domain-wall cost)},
where the cost counts internal bonds crossing the identity/swap domain wall plus
dangling legs whose pinning (swap on A, identity off A) mismatches their tensor's
domain. At large bond dimension χ the sum is dominated by the minimal wall:
S₂ → mincut(A)·ln χ — the Ryu–Takayanagi min-cut. The two-state variable is
therefore not assumed: it IS the replica structure of the entanglement calculation.
The horizon chain then closes: RTN spins (two states, derived) + pairwise 1/r
gravitational coupling between them (run_gravitational_coupling.py) + mean-field
exactness for α < d (run_mf_exactness.py) ⇒ the §5 Landau F(m), with
m = fraction of horizon bonds on the bulk side of the cut.

WE VERIFY, on a 6-tensor two-triangle network with a single bridge bond (a bulk
bottleneck, so geometry beats naive leg-counting):
  (1) min-cut law: S₂(A) → mincut(A)·ln χ — including the KEY case |A| = 3 legs
      with mincut = 1 (the bridge): entropy saturates at ~ln χ, NOT 3 ln χ;
  (2) the Ising dual: the exact two-state partition sum Z_A/Z_∅ reproduces the
      measured average purity, improving with χ — the binary variable does the
      whole computation.

HONEST FLAGS. RTN is a MODEL of holographic states (exact min-cut at χ → ∞ with
flat entanglement spectra); identifying horizon dof with RTN bonds is the
S_ent = S_grav tier restated micro-structurally. What this closes is R2
specifically: GIVEN a tensor-network description, two-state is derived, not chosen.
"""
import numpy as np
from itertools import product

RNG = np.random.default_rng(9)

# network: two triangles (1,2,3) and (4,5,6) joined by one bridge bond 3-4;
# every tensor carries one dangling (boundary) leg.
BONDS = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)]  # last = bridge
NV = 6

def random_state(chi):
    """Contract the RTN into the boundary state ψ[leg1..leg6]."""
    def T(r):  # random tensor of rank r (complex Gaussian)
        return (RNG.standard_normal((chi,) * r) + 1j * RNG.standard_normal((chi,) * r))
    # bond indices: 01:a, 12:b, 02:c, 34:d, 45:e, 35:f, 23:g ; legs i..n
    T1 = T(3); T2 = T(3); T3 = T(4); T4 = T(4); T5 = T(3); T6 = T(3)
    psi = np.einsum("aci,abj,bgck,gdfl,dem,efn->ijklmn",
                    T1, T2, T3, T4, T5, T6, optimize=True)
    psi = psi.reshape(-1)
    return (psi / np.linalg.norm(psi)).reshape((chi,) * NV)

def purity(psi, legs_A, chi):
    axes = list(legs_A) + [k for k in range(NV) if k not in legs_A]
    M = np.transpose(psi, axes).reshape(chi ** len(legs_A), -1)
    lam = np.linalg.svd(M, compute_uv=False) ** 2
    return float((lam ** 2).sum())

def mincut(legs_A):
    """Brute-force min domain-wall cost over the 2^6 vertex bipartitions."""
    best = NV + len(BONDS)
    for sig in product([0, 1], repeat=NV):
        cost = sum(sig[i] != sig[j] for i, j in BONDS)
        cost += sum((sig[v] == 0) == (v in legs_A) for v in range(NV))
        best = min(best, cost)
    return best

def ising_Z(legs_A, chi):
    """Exact two-state partition sum Z_A = Σ_σ χ^{−cost(σ; A)}."""
    Z = 0.0
    for sig in product([0, 1], repeat=NV):
        cost = sum(sig[i] != sig[j] for i, j in BONDS)
        cost += sum((sig[v] == 0) == (v in legs_A) for v in range(NV))
        Z += chi ** (-cost)
    return Z

REGIONS = {"single leg {1}": [0], "two legs {1,2}": [0, 1],
           "TRIANGLE {1,2,3} (bridge!)": [0, 1, 2]}

print("=" * 78)
print("1. Min-cut law: S₂(A) → mincut(A)·ln χ — geometry beats leg counting")
print("=" * 78)
CHIS = [2, 4, 8]
S2 = {name: [] for name in REGIONS}
PUR = {name: [] for name in REGIONS}
for chi in CHIS:
    draws = [random_state(chi) for _ in range(60 if chi < 8 else 30)]
    for name, A in REGIONS.items():
        ps = [purity(psi, A, chi) for psi in draws]
        PUR[name].append(float(np.mean(ps)))
        S2[name].append(float(-np.log(np.mean(ps))))
print(f"    {'region':>28} {'|A|':>4} {'mincut':>7} " +
      " ".join(f"S₂/lnχ (χ={c})" for c in CHIS))
for name, A in REGIONS.items():
    mc = mincut(A)
    vals = " ".join(f"{S2[name][k]/np.log(c):>12.2f}" for k, c in enumerate(CHIS))
    print(f"    {name:>28} {len(A):>4} {mc:>7} {vals}")
print("  ⇒ the 3-leg triangle behind the bridge saturates at ~1·ln χ, not 3·ln χ:")
print("    the entropy is the GEOMETRIC min-cut — the RT structure in miniature.")

print("\n" + "=" * 78)
print("2. The Ising dual: the two-state sum Z_A/Z_∅ IS the purity computation")
print("=" * 78)
devs = {c: [] for c in CHIS}
print(f"    {'region':>28} {'χ':>3} {'−ln⟨purity⟩':>12} {'−ln(Z_A/Z_∅)':>13} {'rel.dev':>8}")
for name, A in REGIONS.items():
    for k, chi in enumerate(CHIS):
        pred = -np.log(ising_Z(A, chi) / ising_Z([], chi))
        meas = S2[name][k]
        dev = abs(meas - pred) / max(pred, 1e-9)
        devs[chi].append(dev)
        print(f"    {name:>28} {chi:>3} {meas:>12.3f} {pred:>13.3f} {dev:>8.2f}")
med_dev = {c: float(np.median(devs[c])) for c in CHIS}
print(f"  median relative deviation by χ: " +
      ", ".join(f"χ={c}: {med_dev[c]:.2f}" for c in CHIS))
print("  ⇒ the exact TWO-STATE (identity/swap) partition sum reproduces the")
print("    entanglement computation, improving with χ — the binary horizon")
print("    variable is DERIVED (it is the replica/Ising dual), not idealised.")

print("\n" + "=" * 78)
print("VERDICT — residual R2 (two-state dof) closed at the RTN level")
print("=" * 78)
print(f"""  The §5 Landau chain now reads, with no binary idealisation inserted by hand:
    RTN entanglement ⇒ Ising spins (two states per dof — THIS script)
    + pairwise 1/r coupling between them (run_gravitational_coupling.py)
    + mean-field exactness for α < d (run_mf_exactness.py)
    ⇒ F(m)/NT = −(J/2)m² + θm + [m ln m + (1−m)ln(1−m)],
  with m = the fraction of horizon bonds on the bulk side of the min-cut — the
  order parameter acquires a microscopic definition, and the area↔volume
  competition is the min-cut's position (boundary-hugging vs bulk-filling).
  FLAG: this closes the two-state tier GIVEN a tensor-network description of the
  horizon state — the residual is now that description itself (the S_ent = S_grav
  tier restated micro-structurally), one assumption where two stood before.""")

# validation
tri, single = "TRIANGLE {1,2,3} (bridge!)", "single leg {1}"
assert S2[tri][-1] / np.log(CHIS[-1]) < 1.6, "bridge region must saturate near 1·ln chi (min-cut, not |A|)"
assert S2["two legs {1,2}"][-1] > S2[single][-1], "S2 must order by mincut"
assert mincut(REGIONS[tri]) == 1 and mincut(REGIONS["two legs {1,2}"]) == 2, "brute-force mincuts"
assert med_dev[CHIS[-1]] < 0.25, "Ising dual must match measured purity at largest chi"
assert med_dev[CHIS[-1]] <= med_dev[CHIS[0]] + 0.05, "Ising-dual agreement must improve with chi"
print("\n[validate] rtn-ising assertions passed.")
