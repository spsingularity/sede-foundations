#!/usr/bin/env python3
"""
Generate the two paper figures that have no existing generator:
  F1 — the f_sat(z) U-shape across cosmic history (Theorem 12), with the fade-in-Λ reading.
  F5 — the false-preference calibration null distribution vs the real Δχ² (Sec. 6.2).
The other four paper figures (F2,F3,F4,F6) reuse results/fig_05, fig_01+fig_03, fig_08, fig_13.
Run: python make_paper_figures.py [n_null]   (default 150)
"""
import sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sede import friedmann as fr

OUT = "output"
Om, gamma = 0.30, 1.4964

def fsat_late(z):
    z = np.atleast_1d(np.asarray(z, float))
    D = fr.compute_growth_factor(z, Om)
    return np.clip((1 - np.exp(-gamma * D**2)) / (1 - np.exp(-gamma)), 0, 1)

# ---------------- F1: f_sat U-shape ----------------
def fig_F1():
    # late-time branch from real growth; early branch is the inflation/reheating de Sitter bracket
    zl = np.logspace(-1, 3.2, 300)          # z = 0.1 .. ~1600 (structure -> recombination)
    fl = fsat_late(zl)
    # epoch anchors for the U (early values are the de Sitter / reheating limits, Thm 2/12)
    anchors_z = np.array([1e27, 1e9, 1100, 1.0, 0.0])
    anchors_f = np.array([0.02, fsat_late(1e9)[0], fsat_late(1100.0)[0], fsat_late(1.0)[0], fsat_late(0.0)[0]])

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    # late branch (rising: structure switches DE on)
    ax.plot(1 + zl, fl, color="C0", lw=2.4, label=r"late branch (structure $\to$ DE), $f_{\rm sat}\uparrow1$")
    # early de Sitter bracket (inflation, f_sat=1) and the reheating sink down to the FRW valley
    ax.plot([1e60, 1e27], [1.0, 1.0], color="C3", lw=2.4, label=r"inflation (de Sitter), $f_{\rm sat}=1$")
    ax.plot([1e27, 1e9], [0.02, fsat_late(1e9)[0]], color="C3", lw=1.6, ls="--",
            label="reheating sink (Thm 2)")
    ax.scatter(1 + anchors_z, anchors_f, color="k", s=16, zorder=5)
    ax.axhline(1.0, color="gray", lw=0.7, ls=":")
    ax.axhline(0.0, color="gray", lw=0.7, ls=":")
    ax.annotate("inflation\n(early bracket)", (3e58, 1.0), xytext=(1e50, 0.80), fontsize=8,
                ha="center", arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.annotate("dark energy\n(late bracket)", (1.0, 1.0), xytext=(8, 0.66), fontsize=8,
                ha="center", arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.annotate("FRW valley\n($f_{\\rm sat}\\approx0$,\nstandard $\\Lambda$CDM era)", (1100, 0.0),
                xytext=(5e4, 0.22), fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.set_xscale("log")
    ax.set_xlabel(r"$1+z$  (cosmic time $\rightarrow$ left)")
    ax.set_ylabel(r"$f_{\rm sat}$  (saturated horizon-entropy fraction)")
    ax.set_title(r"F1 — the $f_{\rm sat}$ U-shape: $\Lambda$ fades in as structure forms")
    ax.set_ylim(-0.05, 1.12)
    ax.invert_xaxis()
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.9)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_F1_fsat_history.png", dpi=130); plt.close(fig)
    print("saved fig_F1_fsat_history.png")

# ---------------- F3: the derived dark-sector inputs ----------------
def fig_F3():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    # left: lambda = 1 - Delta/2
    D = np.linspace(0, 1, 100)
    axL.plot(D, 1 - D / 2, color="C0", lw=2.4)
    axL.scatter([0, 1], [1.0, 0.5], color=["gray", "C3"], s=70, zorder=5)
    axL.annotate(r"$\Delta=0:\ \lambda=1$" + "\n(area-law, disfavoured)", (0, 1.0),
                 xytext=(0.12, 0.86), fontsize=8.5)
    axL.annotate(r"$\Delta=1:\ \lambda=0.5$" + "\n(volume-law, adopted)", (1, 0.5),
                 xytext=(0.30, 0.55), fontsize=8.5, color="C3")
    axL.set_xlabel(r"Barrow deformation $\Delta$")
    axL.set_ylabel(r"H-coupling $\lambda=1-\Delta/2$")
    axL.set_title(r"$\lambda$ from the horizon (Thm 8–9)")
    axL.set_xlim(-0.03, 1.03); axL.set_ylim(0.45, 1.05)
    # right: gamma vs entropy-weight exponent p (binding energy into the horizon)
    from sede.gamma_computation import entropy_weight_scan
    scan = entropy_weight_scan(p_list=(2/3, 1.0, 4/3, 5/3, 2.0))
    ps = np.array([p for p, _ in scan]); gammas = np.array([g for _, g in scan])
    axR.plot(ps, gammas, "o-", color="C0", lw=2)
    # mark the derived p=5/3 -> gamma~1.5
    i53 = int(np.argmin(abs(ps - 5/3)))
    axR.scatter([5/3], [gammas[i53]], color="C3", s=80, zorder=6)
    axR.annotate(fr"$p=5/3 \to \gamma={gammas[i53]:.2f}$" + "\n(binding energy into\nthe horizon, Thm 4C)",
                 (5/3, gammas[i53]), xytext=(0.78, gammas[i53] + 0.05), fontsize=8.5, color="C3")
    axR.set_xlabel(r"entropy-weight exponent $p\quad(\Sigma_S\propto\int M^p\,dn/d\ln M)$")
    axR.set_ylabel(r"$\gamma = d\ln\Sigma_S/d\ln\sigma_8$")
    axR.set_title(r"$\gamma$ from structure (Thm 4C)")
    fig.suptitle("F3 — the derived dark-sector inputs: no fitted parameter", y=1.02)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_F3_derived_inputs.png", dpi=130,
                                    bbox_inches="tight"); plt.close(fig)
    print(f"saved fig_F3_derived_inputs.png  (gamma(p=5/3)={gammas[i53]:.3f})")

# ---------------- F5: false-preference calibration ----------------
def fig_F5(n_null=150):
    from run_xval_calibration import make_mock, fit_mock, data, CMB_CINV
    from scipy.linalg import cholesky
    import sede.data_loader as dl
    _, _, covp = dl.load_pantheon_plus(); data["_panL"] = cholesky(covp, lower=True)
    data["_cmbL"] = cholesky(np.linalg.inv(CMB_CINV), lower=True)
    truth_L = (0.299, 68.66, 0.02237, -19.40, 0.760)
    REAL = -2.96
    d = []
    for i in range(n_null):
        mk = make_mock(truth_L, True)
        d.append(fit_mock(False, mk) - fit_mock(True, mk))
        if (i + 1) % 25 == 0: print(f"  [null] {i+1}/{n_null} ...", flush=True)
    d = np.array(d)
    k = int(np.sum(d <= REAL)); p = k / n_null
    np.save(f"{OUT}/calib_null.npy", d)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.hist(d, bins=24, color="C0", alpha=0.8, edgecolor="white",
            label=fr"$\Lambda$CDM-truth null ($N={n_null}$), mean$={d.mean():+.2f}$")
    ax.axvline(REAL, color="C3", lw=2.4, label=fr"real data $\Delta\chi^2={REAL:+.2f}$")
    ax.axvline(0, color="gray", lw=0.8, ls=":")
    ax.set_xlabel(r"$\Delta\chi^2(\mathrm{SEDE}-\Lambda\mathrm{CDM})$")
    ax.set_ylabel("mocks")
    ax.set_title("F5 — false-preference calibration: the preference is not flexibility")
    ax.text(0.02, 0.97, f"real beyond {k}/{n_null} of the null\n"
                        f"(headline 500-mock run: $p\\approx0.002$, ~3$\\sigma$)\n"
                        f"this run: $p={p:.3f}$",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_F5_calibration.png", dpi=130); plt.close(fig)
    print(f"saved fig_F5_calibration.png  (null mean {d.mean():+.2f}, p={p:.3f})")

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    fig_F1()
    fig_F5(n)
