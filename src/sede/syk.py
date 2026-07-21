"""
SYK as the horizon scrambler — a literature-grade scrambling experiment for SEDE's
volume-law (Δ=1) postulate.

Context. SEDE's open postulate (Thm 10) is that the cosmic horizon is in a maximally
*thermalised / scrambled* state, whose entanglement is therefore VOLUME-law (S ∝ V ∝
A^{3/2} ⟹ Δ=1) rather than the ground-state AREA-law. The sibling team's S6 toy showed
the area→volume transition with a generic random-circuit drive (run_scrambling_toy.py).
This module replaces that ad-hoc circuit with the **Sachdev–Ye–Kitaev model** — the
canonical maximally-chaotic, *holographically dual* (near-extremal black hole) model
that *saturates* the Maldacena–Shenker–Stanford chaos bound λ_L ≤ 2πT. We then run three
*respected* diagnostics, each a standard certificate in the quantum-chaos literature:

  (1) Spectral statistics — the gap ratio ⟨r⟩ reproduces Wigner–Dyson RMT, including the
      famous N mod 8 symmetry class (GOE/GUE/GSE; García-García & Verbaarschot 2016;
      Cotler et al 2017). This is the operational meaning of "maximal scrambler".
  (2) Eigenstate entanglement — mid-spectrum SYK eigenstates saturate the Page (volume-
      law) value (Vidmar–Rigol). This is the Δ=0→Δ=1 (area→volume) content SEDE needs,
      in genuine eigenstates of the black-hole-dual model — not a hand-built circuit.
  (3) OTOC dynamics — the infinite-temperature out-of-time-order correlator scrambles
      completely (C(t)→2, the maximal value) on a short time, vs an integrable control.

Control in every diagnostic: the q=2 (free-Majorana, quadratic) model, which is
integrable — Poisson-trending statistics, sub-maximal entanglement, incomplete scrambling.

Honest scope: this is a many-body quantum model, not gravity. It does NOT derive Δ=1 for
the real horizon or fix the dark-energy magnitude. What it establishes rigorously is the
*physics SEDE invokes*: a maximally-scrambled (SYK-class) system has volume-law (Δ=1)
entanglement. The one open step stays open — that the COSMIC horizon is held in that
driven steady state (Route 2 structure deposition; the MSS bound gives the scrambling
time (ln S)/H ≈ 282/H reproduced in scrambling_time_bridge()).
"""
import math
import numpy as np
from itertools import combinations

# ── Pauli matrices ────────────────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Gap-ratio references (Atas et al 2013) and Poisson; SYK class by N mod 8.
R_REF = {"Poisson": 0.3863, "GOE": 0.5307, "GUE": 0.5996, "GSE": 0.6744}
SYK_CLASS = {0: "GOE", 2: "GUE", 4: "GSE", 6: "GUE"}


def _kron(ops):
    out = np.array([[1]], dtype=complex)
    for o in ops:
        out = np.kron(out, o)
    return out


def majoranas(N):
    """N Majorana operators (N even) on n=N/2 qubits via Jordan–Wigner.
    Normalisation {γ_a, γ_b} = 2δ_ab (γ_a²=I, Hermitian):
        γ_{2k}   = Z…Z X,  γ_{2k+1} = Z…Z Y  (Z on qubits 0..k-1)."""
    assert N % 2 == 0, "N must be even"
    n = N // 2
    g = []
    for k in range(n):
        zstring = [Z] * k
        for P in (X, Y):
            g.append(_kron(zstring + [P] + [I2] * (n - k - 1)))
    return g


def syk_hamiltonian(N, q=4, rng=None, J=1.0):
    """Majorana SYK_q Hamiltonian (q=4 chaotic / q=2 free), real Gaussian couplings.
       q=4: H = Σ_{i<j<k<l} J_{ijkl} γ_iγ_jγ_kγ_l,  Var = 3! J²/N³
       q=2: H = Σ_{i<j} i K_{ij} γ_iγ_j,            Var = J²/N (free/integrable)."""
    rng = np.random.default_rng() if rng is None else rng
    g = majoranas(N)
    D = g[0].shape[0]
    H = np.zeros((D, D), dtype=complex)
    if q == 4:
        sig = np.sqrt(math.factorial(3) * J**2 / N**3)
        for (i, j, k, l) in combinations(range(N), 4):
            H += rng.normal(0.0, sig) * (g[i] @ g[j] @ g[k] @ g[l])
    elif q == 2:
        sig = np.sqrt(J**2 / N)
        for (i, j) in combinations(range(N), 2):
            H += rng.normal(0.0, sig) * 1j * (g[i] @ g[j])
    else:
        raise ValueError(f"q={q} not supported (use 2 or 4)")
    return 0.5 * (H + H.conj().T), g


def parity_sector_indices(N, sector=+1):
    """Computational-basis indices with given fermion parity P=∏Z (H is block-diagonal
    in parity since q is even). Restricting avoids mixing symmetry sectors."""
    n = N // 2
    return np.array([i for i in range(2**n)
                     if (bin(i).count("1") % 2 == 0) == (sector == +1)])


# ── (1) spectral statistics ─────────────────────────────────────────────────────
def gap_ratio(evals, collapse_tol=1e-9):
    """Mean gap ratio ⟨r_n⟩, r_n=min(s_n,s_{n+1})/max(...), central 80%. Collapses exact
    degeneracies first (needed for GSE Kramers doubling)."""
    e = np.sort(evals.real)
    e = e[np.concatenate([[True], np.diff(e) > collapse_tol])]
    s = np.diff(e)
    s = s[s > 0]
    if len(s) < 3:
        return np.nan
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    lo, hi = int(0.1 * len(r)), int(0.9 * len(r))
    return float(np.mean(r[lo:hi]))


def mean_gap_ratio(N, q, n_real, seed=0):
    """Disorder-averaged ⟨r⟩ in the even-parity sector (mean, standard error)."""
    rng = np.random.default_rng(seed)
    idx = parity_sector_indices(N, +1)
    rs = []
    for _ in range(n_real):
        H, _ = syk_hamiltonian(N, q=q, rng=rng)
        r = gap_ratio(np.linalg.eigvalsh(H[np.ix_(idx, idx)]))
        if not np.isnan(r):
            rs.append(r)
    return float(np.mean(rs)), float(np.std(rs) / np.sqrt(max(len(rs), 1)))


# ── (2) eigenstate entanglement (Page / volume law) ─────────────────────────────
def entanglement_entropy(psi, n, n_A):
    """von Neumann entropy (nats) of the first n_A qubits of state vector psi (dim 2^n)."""
    M = psi.reshape(2**n_A, 2**(n - n_A))
    p = np.linalg.svd(M, compute_uv=False) ** 2
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def page_value(n_A, n):
    """Haar-random pure-state average entanglement (nats), Page's formula:
    ⟨S⟩ = Σ_{k=d_B+1}^{d_A d_B} 1/k − (d_A−1)/(2 d_B),  d_A=2^{n_A} ≤ d_B."""
    dA, dB = 2**n_A, 2**(n - n_A)
    if dA > dB:
        dA, dB = dB, dA
    return float(np.sum(1.0 / np.arange(dB + 1, dA * dB + 1)) - (dA - 1.0) / (2.0 * dB))


def page_curve(N, q, n_real, n_states=6, seed=1):
    """Mean mid-spectrum eigenstate entanglement S(n_A) vs Page; returns (n_A, S, S_page)."""
    rng = np.random.default_rng(seed)
    n = N // 2
    nAs = list(range(1, n))
    acc = {a: [] for a in nAs}
    for _ in range(n_real):
        H, _ = syk_hamiltonian(N, q=q, rng=rng)
        ev, U = np.linalg.eigh(H)
        mid = len(ev) // 2
        for s in range(mid - n_states // 2, mid + (n_states + 1) // 2):
            for a in nAs:
                acc[a].append(entanglement_entropy(U[:, s], n, a))
    return (np.array(nAs),
            np.array([np.mean(acc[a]) for a in nAs]),
            np.array([page_value(a, n) for a in nAs]))


# ── (3) OTOC scrambling dynamics ────────────────────────────────────────────────
def otoc_curve(N, q, ts, n_real, seed=2):
    """Infinite-T OTOC F(t)=⟨W(t)V W(t)V⟩/D, C(t)=2(1−Re F), with Hermitian unitary
    DISJOINT bilinears W=iγ_0γ_1, V=iγ_{N-2}γ_{N-1} ([W,V]=0 at t=0 ⟹ C(0)=0)."""
    rng = np.random.default_rng(seed)
    Fs = np.zeros((n_real, len(ts)))
    for r in range(n_real):
        H, g = syk_hamiltonian(N, q=q, rng=rng)
        D = H.shape[0]
        W = 1j * g[0] @ g[1]
        V = 1j * g[N - 2] @ g[N - 1]
        E, U = np.linalg.eigh(H)
        Wd, Vd = U.conj().T @ W @ U, U.conj().T @ V @ U
        for it, t in enumerate(ts):
            ph = np.exp(1j * E * t)
            Wt = (ph[:, None] * Wd) * ph.conj()[None, :]
            Fs[r, it] = np.trace(Wt @ Vd @ Wt @ Vd).real / D
    F = Fs.mean(axis=0)
    return F, 2.0 * (1.0 - F)


# ── (4) the MSS chaos-bound → SEDE scrambling-time bridge ────────────────────────
def scrambling_time_bridge(lnS=281.73):
    """SYK saturates the MSS bound λ_L = 2πT (fastest scrambler). A scrambler reaches
    volume-law in t* = (ln S)/λ_L = (ln S)/(2πT). For the cosmic horizon T_dS=H/2π ⟹
    2πT = H ⟹ **t* = (ln S)/H**. With horizon ln S ≈ 282 (Route 2) this is the ~282/H
    scrambling time that pins the driven volume-law steady state. Returns the dict."""
    return {"lambda_L_over_2piT": 1.0,            # MSS saturation (SYK)
            "lnS": lnS,
            "t_star_over_inv_H": lnS,             # t* = (ln S)/H
            "route2_lnS": 281.73}
