"""
Tensor-network (Ryu–Takayanagi min-cut) test of what makes a horizon area- vs volume-law.

For a random tensor network the entanglement of a region is S = (min-cut separating it from the
rest) × ln χ (Ryu–Takayanagi; χ = bond dim, take ln χ = ln 2). We place N nodes on a D-dim spatial
grid, enclose a ball of radius R, and compute the min-cut (= max-flow) between the ball interior and
the exterior under two connectivities:

  LOCAL    (nearest-neighbour lattice): min-cut ∝ boundary ∝ R^{D-1}  -> AREA law (Δ=0)
  NONLOCAL (random-regular / expander): min-cut ∝ interior size ∝ R^D -> VOLUME law (Δ=1)

So volume-law (Δ=1) horizon entropy is REALISABLE in a concrete holographic-code structure — its
required input is NONLOCAL connectivity of the horizon dof (the geometry-free, all-to-all structure
that makes SYK maximally chaotic). A geometrically local horizon is area-law. This is an existence
proof for the counting half of the postulate, and names the needed ingredient: nonlocality.
(2D testbed: AREA ∝ R¹, VOLUME ∝ R²; in physical 3D space these are R² vs R³ = Δ=0 vs Δ=1.)
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

BIG = 10**7


def _grid(n):
    """N=n² node positions on an n×n grid; returns positions and index lookup."""
    xs, ys = np.meshgrid(np.arange(n), np.arange(n))
    pos = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    return pos


def _local_edges(n):
    """nearest-neighbour (4-neighbour) lattice edges as (i,j) pairs."""
    edges = []
    idx = lambda a, b: a * n + b
    for a in range(n):
        for b in range(n):
            if a + 1 < n:
                edges.append((idx(a, b), idx(a + 1, b)))
            if b + 1 < n:
                edges.append((idx(a, b), idx(a, b + 1)))
    return edges


def _nonlocal_edges(N, deg, rng):
    """random deg-regular-ish graph: each node linked to `deg` random others (position-blind)."""
    edges = set()
    for i in range(N):
        for _ in range(deg):
            j = int(rng.integers(N))
            if j != i:
                edges.add((min(i, j), max(i, j)))
    return list(edges)


def _min_cut(N, edges, interior_mask):
    """min-cut (= max-flow) separating interior nodes from exterior, unit-capacity edges.
    Super-source -> interior (BIG), exterior -> super-sink (BIG)."""
    src, snk = N, N + 1
    M = N + 2
    rows, cols, dat = [], [], []
    for (i, j) in edges:                       # undirected unit edges
        rows += [i, j]; cols += [j, i]; dat += [1, 1]
    for i in np.where(interior_mask)[0]:       # source -> interior
        rows.append(src); cols.append(i); dat.append(BIG)
    for i in np.where(~interior_mask)[0]:      # exterior -> sink
        rows.append(i); cols.append(snk); dat.append(BIG)
    g = csr_matrix((np.array(dat), (np.array(rows), np.array(cols))), shape=(M, M))
    return int(maximum_flow(g, src, snk).flow_value)


def mincut_scaling(n=44, deg=6, Rs=(4, 6, 8, 10, 12), seed=0):
    """min-cut S(R) for LOCAL vs NONLOCAL connectivity on the same N=n² nodes."""
    rng = np.random.default_rng(seed)
    pos = _grid(n)
    N = n * n
    c = np.array([(n - 1) / 2.0, (n - 1) / 2.0])
    rad = np.hypot(pos[:, 0] - c[0], pos[:, 1] - c[1])
    le = _local_edges(n)
    ne = _nonlocal_edges(N, deg, rng)
    Rs = list(Rs)
    loc = [_min_cut(N, le, rad < R) for R in Rs]
    non = [_min_cut(N, ne, rad < R) for R in Rs]
    fit = lambda y: float(np.polyfit(np.log(Rs), np.log(np.maximum(y, 1)), 1)[0])
    return dict(n=n, N=N, deg=deg, Rs=Rs, local=loc, nonlocal_=non,
                exp_local=fit(loc), exp_nonlocal=fit(non))
