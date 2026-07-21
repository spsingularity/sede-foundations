"""
make_foundations_figures.py — figures for the SEDE foundations paper (SEDE_foundations.md).

Produces the manuscript's four figures into results/:
  foundations_fig1_reduction.png   — Fig 1, the reformulation-chain schematic (the spine)
  foundations_fig_birth.png        — Fig 4, birth at the ordering bifurcation + barrier scaling
  foundations_fig5_inferences.png  — Fig 2, ensemble inequivalence (C<0) + GREA viscosity
  criticality_soc.png              — Fig 3, produced by experiments/run_chr_soc.py ('soc' stage)

The SOC-signature figure (foundations_fig4_socsignature.png) is retired: the SOC/roughening
z* layer it depicted was retracted (§6).

Reproduced from the single entry point:  python reproduce.py --only foundfig,soc
(or directly: python figures/make_foundations_figures.py && python experiments/run_chr_soc.py)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "output")
os.makedirs(OUT, exist_ok=True)

GREEN, BLUE, ORANGE, GREY = "#2e7d32", "#1565c0", "#e65100", "#555555"


# ---------------------------------------------------------------------------
# Fig 1 — the reduction chain (schematic)
# ---------------------------------------------------------------------------
def fig_reduction():
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    def box(x, y, w, h, text, color, fc=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                     linewidth=1.6, edgecolor=color, facecolor=fc or "white"))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=8.6, color="black")

    def arrow(x1, y1, x2, y2, color=GREY):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=13, linewidth=1.3, color=color))

    box(3.3, 9.0, 3.4, 0.8, "Volume-law postulate\n(horizon counts its volume, Δ=1)", GREY, "#f2f2f2")
    # four sub-claims (state/scale are conditional/reduced, not outright derived — see §3)
    subs = [("state\nSYK (route A, conditional)", GREEN, 0.3),
            ("form\nreduced (thermal)", BLUE, 2.7),
            ("scale\nreduced (CKN bound)", BLUE, 5.1),
            ("count\nOPEN — the residue", ORANGE, 7.5)]
    for txt, col, x in subs:
        box(x, 7.4, 2.2, 0.85, txt, col)
        arrow(5.0, 9.0, x + 1.1, 8.25, col)
    # chain from count
    arrow(8.6, 7.4, 8.6, 6.7, ORANGE)
    box(5.6, 5.9, 4.0, 0.8, "count ⟺ connectivity ⟸ range:\ngravity 1/r (α≤d) ⇒ non-additive ⇒ volume class", BLUE)
    arrow(7.6, 5.9, 6.2, 5.3, BLUE)
    box(3.0, 4.5, 5.0, 0.8, "F(m) derived given the two-state dof idealisation:\nexact counting + pairwise-gravity m² + mean-field (α≤d)", GREEN)
    arrow(5.5, 4.5, 5.5, 3.9, GREEN)
    box(2.8, 3.1, 5.4, 0.8, "BIRTH SELECTION: J∝N crosses J_c at N~1; Barrow-tilted\nbifurcation ⇒ volume (P≥0.9); hysteresis ⇒ BH area", GREEN)
    arrow(5.5, 3.1, 5.5, 2.5, GREY)
    box(3.2, 1.7, 4.6, 0.8, "RESIDUE: FRW-patch algebra type (II₁ vs II_∞) —\nits trace calibration is the Δ measurement (§7.1b)", ORANGE)
    arrow(5.5, 1.7, 5.5, 1.1, GREEN)
    box(2.6, 0.3, 5.8, 0.7, "DECIDED EMPIRICALLY: point prediction Δ = 1 (ceiling Δ ≤ 1) — DESI DR3 + Euclid (~0.09)", GREEN, "#eaf5ea")

    # legend
    for i, (c, t) in enumerate([(GREEN, "derived / decided"), (BLUE, "reduced"), (ORANGE, "open")]):
        ax.add_patch(FancyBboxPatch((0.2, 0.3 + i*0.5), 0.3, 0.3, boxstyle="round,pad=0.02",
                     edgecolor=c, facecolor="white", linewidth=1.6))
        ax.text(0.6, 0.45 + i*0.5, t, fontsize=7.5, va="center")
    ax.set_title("The reformulation of SEDE's one postulate (birth-selection chain)", fontsize=11)
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig1_reduction.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------------------
# Fig birth — selection at the ordering bifurcation (run_birth_bifurcation.py)
# ---------------------------------------------------------------------------
def fig_birth():
    B = 0.3                                              # Barrow tilt per dof
    U = np.linspace(-25, 25, 40001); M = 1.0/(1.0 + np.exp(-U))
    def f(m, J, th): return -0.5*J*m**2 + th*m + m*np.log(m) + (1-m)*np.log(1-m)
    def fp(m, J, th): return -J*m + th + np.log(m/(1-m))
    def minima(J, th):
        y = f(M, J, th)
        i = np.where((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:]))[0] + 1
        return M[i]

    # (a) tilted-bifurcation branch diagram + swept gradient-flow trajectory
    Js = np.linspace(2.0, 24.0, 400); vol, area = [], []
    for J in Js:
        mn = minima(J, J/2 - B)
        vol.append(mn.max()); area.append(mn.min() if len(mn) > 1 else np.nan)
    vol, area = np.array(vol), np.array(area)
    dt, steps, cap = 0.005, 100000, 2e-3
    Jp = np.linspace(2.0, 24.0, steps); m = 0.45; traj = []
    for i in range(steps):
        J = Jp[i]
        dm = float(np.clip(-fp(m, J, J/2 - B)*dt, -cap, cap))
        m = float(np.clip(m + dm, 1e-9, 1 - 1e-9))
        if i % 500 == 0: traj.append((J, m))
    traj = np.array(traj)

    # (b) barrier vanishes quadratically at threshold (coexistence tilt b=0)
    dJ = np.array([0.1, 0.2, 0.4, 0.8, 1.2]); bars = []
    for d in dJ:
        J = 4.0 + d; y = f(M, J, J/2); mn = minima(J, J/2)
        lo, hi = mn.min(), mn.max()
        bars.append(y[(M > lo) & (M < hi)].max() - f(lo, J, J/2))
    bars = np.array(bars)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(Js, vol, color=GREEN, lw=2.2, label="volume branch (connected to disorder)")
    ax[0].plot(Js, area, color=ORANGE, lw=2, ls="--", label="area branch (born metastable)")
    ax[0].plot(traj[:, 0], traj[:, 1], color=BLUE, lw=1.1, alpha=0.9,
               label="grown horizon (gradient flow)")
    ax[0].axvline(4.0, ls=":", c="0.6"); ax[0].text(4.1, 0.03, r"$J_c$", fontsize=9, color=GREY)
    J_sp = float(Js[np.isfinite(area)].min())
    ax[0].axvline(J_sp, ls=":", c=ORANGE, alpha=0.7)
    ax[0].text(J_sp + 0.1, 0.16, r"$J_{sp}$", fontsize=9, color=ORANGE)
    ax[0].annotate("BH: created abruptly at large N,\nlocked by $e^{-N\\Delta f}$ (hysteresis)",
                   xy=(20, float(area[np.isfinite(area)][-1]) + 0.02), xytext=(11.5, 0.32),
                   fontsize=7.5, color=GREY, arrowprops=dict(arrowstyle="->", color="0.6"))
    ax[0].set_xlabel(r"$J/T \propto 2\pi G m^2 N$  (grows with the horizon)")
    ax[0].set_ylabel("order parameter m (volume-participating fraction)")
    ax[0].set_title("(a) selection at the Barrow-tilted ordering bifurcation")
    ax[0].legend(fontsize=7.5, loc="center right")
    ax[1].loglog(dJ, bars, "o-", color=BLUE, label="barrier $\\Delta f$ (coexistence)")
    ax[1].loglog(dJ, (3/64)*dJ**2, "k--", label=r"Landau: $\frac{3}{64}(J-J_c)^2$")
    ax[1].set_xlabel(r"$J - J_c$"); ax[1].set_ylabel(r"barrier per dof $\Delta f/T$")
    ax[1].set_title("(b) barrier vanishes at threshold — nothing to cross")
    ax[1].legend(fontsize=8)
    ax[1].text(0.05, 0.62, "today: $N\\,\\Delta f \\sim 10^{122}$\n⇒ birth branch locked forever",
               transform=ax[1].transAxes, fontsize=7.5, color=GREY)
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig_birth.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------------------
# Fig 4 — mean-field SOC signature
# ---------------------------------------------------------------------------
def fig_socsignature():
    rng = np.random.default_rng(0)
    def avalanche(rng, eps=0.0, cap=200000):
        active, size = 1, 0
        while active > 0 and size < cap:
            b = rng.poisson((1.0 - eps) * active); size += active; active = b
        return size
    def pdf_of(eps, n=80000):
        sizes = np.array([avalanche(rng, eps) for _ in range(n)])
        sizes = sizes[(sizes > 0) & (sizes < 200000)]
        edges = np.unique(np.round(np.logspace(0, np.log10(sizes.max()), 24)).astype(int))
        cnt, _ = np.histogram(sizes, bins=edges)
        ctr = np.sqrt(edges[:-1]*edges[1:]); pdf = cnt/np.diff(edges)/cnt.sum()
        ok = pdf > 0
        return ctr[ok], pdf[ok]
    # derived dissipation ε = |κ|·t_dyn with exact κ = H(1−3w_eff)/4 at z≈0 (w_eff≈−0.7)
    EPS_H = (1 - 3 * (-0.7)) / 4.0 * 178.0 ** -0.5
    ctr, pdf = pdf_of(0.0)                               # conservative reference
    ctr_h, pdf_h = pdf_of(EPS_H)                         # horizon (truncated)
    ok = np.ones(len(pdf), bool)

    from sede import chr_mechanism as chr
    from scipy.integrate import solve_ivp
    J, TH = 6.0, 3.0; sig = lambda x: 1/(1+np.exp(-x))
    z_star = chr.transition_redshift()
    ng = np.linspace(-4, 1.5, 700); zg = np.exp(-ng)-1
    fz = np.clip(chr.f_eq(chr.control_variance(np.clip(zg, 0, None))), 1e-6, 1)
    def rhs(n, m):
        m = float(np.clip(np.ravel(m)[0], 0, 1))   # solve_ivp passes a 1-elem array (numpy>=2.5-safe)
        return 25*(-m+sig(J*m+float(np.interp(n, ng, fz))-TH))
    m = solve_ivp(rhs, (ng[0], ng[-1]), [sig(-TH)], t_eval=ng, rtol=1e-8, atol=1e-10).y[0]
    arg = J*m+fz-TH; lam = 1-J*sig(arg)*(1-sig(arg)); chi = 1/np.maximum(np.abs(lam), 1e-3)
    m_sp = 0.5 - np.sqrt(0.25 - 1.0/J)               # area-branch spinodal
    chi_area = np.where(m <= m_sp, chi, np.nan)       # show only the metastable area branch

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(ctr, pdf, "o", ms=4, color=BLUE, label=r"conservative ($\epsilon\to 0$)")
    ax[0].loglog(ctr_h, pdf_h, "s", ms=4, color=ORANGE,
                 label=rf"horizon $\epsilon = |\kappa|\,t_{{dyn}} \approx {EPS_H:.3f}$ (z=0; smaller at z*)")
    s = np.array([2, ctr.max()/15]); ax[0].loglog(s, 0.5*s**-1.5, "k--", label=r"$P(s)\propto s^{-3/2}$")
    s_c = 2.0 / EPS_H**2
    sfull = np.geomspace(2, 5*s_c, 100)
    ax[0].loglog(sfull, 0.5*sfull**-1.5*np.exp(-sfull/s_c), ":", color=ORANGE, lw=1.5,
                 label=rf"$s^{{-3/2}}e^{{-s/s_c}}$, $s_c=2/\epsilon^2\approx{s_c:.0f}$")
    ax[0].axvline(s_c, ls=":", c="0.75", lw=1)
    ax[0].set_xlabel("avalanche size s"); ax[0].set_ylabel("P(s)")
    ax[0].set_title("(a) SOC signature — truncated window (membrane-derived)")
    ax[0].legend(fontsize=7.5)
    ax[0].text(0.04, 0.06, "nearly conservative automatically:\nfast (t_dyn) redistribution, Hubble-slow sink",
               transform=ax[0].transAxes, fontsize=7, color=GREY, va="bottom")
    sel = (zg > 0) & (zg < 3.5)
    ax[1].plot(zg[sel], chi_area[sel], color=ORANGE, lw=2)
    ax[1].axvline(z_star, ls=":", c="0.6", label=f"z*={z_star:.2f}")
    ax[1].set_xlim(3.5, 0); ax[1].set_xlabel("redshift z")
    ax[1].set_ylabel(r"susceptibility / $\tau_{ac}$  ($1/\lambda_{stab}$)")
    ax[1].set_title("(b) transient critical slowing at z*"); ax[1].legend(fontsize=8)
    ax[1].annotate("area-phase response\npeaks with deposition", xy=(z_star, 0.85*np.nanmax(chi_area[sel])),
                   xytext=(2.3, 0.6*np.nanmax(chi_area[sel])), fontsize=7.5, color=GREY,
                   arrowprops=dict(arrowstyle="->", color="0.6"))
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig4_socsignature.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------------------
# Fig 5 — literature inferences: ensemble (C<0) + GREA viscosity
# ---------------------------------------------------------------------------
def fig_inferences():
    from sede.friedmann import E_SEDE_lambda, compute_growth_factor
    Om, Or, g, lam = 0.30, 9e-5, 1.4964, 0.5; ODE0 = 1-Om-Or
    z = np.linspace(0, 3, 300); H = E_SEDE_lambda(z, Om, g, lam, Or)
    S = 1/H**2; T = H/(2*np.pi); C = T*np.gradient(S, T)
    a = 1/(1+z); D = compute_growth_factor(z, Om, Or)
    f = np.clip((1-np.exp(-g*D**2))/(1-np.exp(-g)), 1e-9, 1)
    lna = np.log(a); eps = -np.gradient(np.log(H), lna); dlnf = np.gradient(np.log(f), lna)
    rho = ODE0*f*H**(2*lam); zeta = rho*dlnf/(9*H**2); w = -1+(1/3)*(2*lam*eps-dlnf)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(z, C/S, color=BLUE, lw=2)
    ax[0].axhline(-2, ls="--", c="0.6", label="C = −2S (analytic)")
    ax[0].set_xlabel("redshift z"); ax[0].set_ylabel("C / S")
    ax[0].set_title("(a) horizon negative specific heat\n⇒ microcanonical (volume) selected")
    ax[0].legend(fontsize=8); ax[0].set_ylim(-2.6, 0)
    ax[0].axhspan(-2.6, 0, color="#fdecea", alpha=0.5, zorder=0)
    ax[0].text(1.5, -0.45, "C < 0  ⇒  canonical ensemble ill-defined\n(ensemble inequivalence, long-range)",
               ha="center", fontsize=7.8, color=GREY)
    ax2 = ax[1]; ln1 = ax2.plot(z, zeta, color=ORANGE, lw=2, label=r"$\zeta(z)$ (bulk viscosity)")
    ax2.set_xlabel("redshift z"); ax2.set_ylabel(r"effective bulk viscosity $\zeta$", color=ORANGE)
    ax2.tick_params(axis="y", labelcolor=ORANGE); ax2.set_xlim(3, 0)
    axb = ax2.twinx(); ln2 = axb.plot(z, w, color=GREY, lw=1.5, ls="--", label="w(z)")
    axb.axhline(-1, ls=":", c="0.7"); axb.set_ylabel("w(z)", color=GREY)
    ax2.set_title("(b) SEDE as structure-gated GREA\n(ζ≥0, peaks at gate-activation z≲0.5)")
    lns = ln1+ln2; ax2.legend(lns, [l.get_label() for l in lns], fontsize=8, loc="center right")
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig5_inferences.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


if __name__ == "__main__":
    fig_reduction()
    fig_birth()
    fig_inferences()
    # fig_socsignature() — retired: the SOC/roughening z* layer was retracted (§6), see D1.
    print("Figs 1, 2, 4 done. For Fig 3 (criticality_soc.png) run the 'soc' stage: python experiments/run_chr_soc.py")
