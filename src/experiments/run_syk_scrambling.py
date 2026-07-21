"""
SYK scrambling experiment — a literature-grade version of the area→volume (Δ=0→Δ=1)
demonstration SEDE's volume-law postulate (Thm 10) requires.

Upgrades the sibling team's S6 random-circuit toy (run_scrambling_toy.py) to the
Sachdev–Ye–Kitaev model: the canonical maximally-chaotic, holographically-dual
(near-extremal black hole) model that saturates the MSS chaos bound λ_L ≤ 2πT. Runs
three respected quantum-chaos diagnostics, each with an integrable q=2 control:

  (1) gap-ratio ⟨r⟩ → Wigner–Dyson RMT incl. the N mod 8 class (GOE/GUE/GSE)
  (2) eigenstate entanglement → Page (volume-law) saturation
  (3) infinite-T OTOC → complete fast scrambling (C→2)

then bridges the MSS bound to SEDE's (ln S)/H ≈ 282/H scrambling time (Route 2).

Output: console summary, results/syk_scrambling.json, results/fig_syk_scrambling.png.
"""
import json
import os
import numpy as np

from sede import syk

# realisation counts (modest defaults; raise for publication-smooth curves)
N_GAP = 250
N_PAGE = 40
N_OTOC = 40


def main():
    out = {"diagnostics": {}}

    # ── (1) gap-ratio RMT statistics ────────────────────────────────────────────
    print("=" * 74)
    print("(1) SPECTRAL STATISTICS — gap ratio ⟨r⟩  (even-parity sector)")
    print("    refs: Poisson=%.3f  GOE=%.3f  GUE=%.3f  GSE=%.3f"
          % (syk.R_REF["Poisson"], syk.R_REF["GOE"], syk.R_REF["GUE"], syk.R_REF["GSE"]))
    print("=" * 74)
    gap = {}
    for N in (8, 10, 12):
        r4, e4 = syk.mean_gap_ratio(N, 4, N_GAP, seed=N)
        r2, e2 = syk.mean_gap_ratio(N, 2, N_GAP, seed=N + 100)
        cls = syk.SYK_CLASS[N % 8]
        closest = min(syk.R_REF, key=lambda k: abs(syk.R_REF[k] - r4))
        ok = closest == cls
        print(f"  N={N:2d} (mod8={N%8}, predict {cls}): q=4 ⟨r⟩={r4:.3f}±{e4:.3f} → {closest} "
              f"{'✓' if ok else '✗'}   | q=2 ⟨r⟩={r2:.3f} (→Poisson)")
        gap[N] = dict(r4=r4, e4=e4, r2=r2, e2=e2, predict=cls, closest=closest, ok=ok)
    out["diagnostics"]["gap_ratio"] = gap

    # ── (2) Page-curve eigenstate entanglement ──────────────────────────────────
    print("\n" + "=" * 74)
    print("(2) EIGENSTATE ENTANGLEMENT — Page (volume-law) saturation")
    print("=" * 74)
    page = {}
    for N in (10, 12):
        nAs, S4, Sp = syk.page_curve(N, 4, N_PAGE, seed=N)
        _, S2, _ = syk.page_curve(N, 2, N_PAGE, seed=N + 7)
        ratio = (S4 / Sp)
        print(f"  N={N} (n={N//2} qubits): q=4/Page (per n_A) = "
              + " ".join(f"{v:.2f}" for v in ratio)
              + f"   [mean {ratio.mean():.3f}]")
        print(f"           q=2 free / Page              = "
              + " ".join(f"{v:.2f}" for v in (S2 / Sp)))
        page[N] = dict(n_A=nAs.tolist(), S_q4=S4.tolist(), S_page=Sp.tolist(),
                       S_q2=S2.tolist(), ratio_q4=ratio.tolist(),
                       mean_ratio_q4=float(ratio.mean()),
                       mean_ratio_q2=float((S2 / Sp).mean()))
    out["diagnostics"]["page_curve"] = page

    # ── (3) OTOC scrambling dynamics ────────────────────────────────────────────
    print("\n" + "=" * 74)
    print("(3) OTOC DYNAMICS — complete fast scrambling (C→2)")
    print("=" * 74)
    ts = np.linspace(0, 8, 41)
    F4, C4 = syk.otoc_curve(12, 4, ts, N_OTOC, seed=12)
    F2, C2 = syk.otoc_curve(12, 2, ts, N_OTOC, seed=62)
    plateau = float(C4[len(C4) // 2:].mean())
    tstar = float(np.interp(0.5 * plateau, C4, ts))
    print(f"  N=12: q=4 C(0)={C4[0]:.2f} → plateau {plateau:.2f} (max=2, full scramble), t*={tstar:.2f}")
    print(f"        q=2 free plateau {float(C2[len(C2)//2:].mean()):.2f} (incomplete scrambling)")
    out["diagnostics"]["otoc"] = dict(ts=ts.tolist(), C_q4=C4.tolist(), C_q2=C2.tolist(),
                                      plateau_q4=plateau, t_star=tstar)

    # ── (4) MSS → SEDE scrambling-time bridge ───────────────────────────────────
    print("\n" + "=" * 74)
    print("(4) MSS CHAOS BOUND → SEDE scrambling time")
    print("=" * 74)
    br = syk.scrambling_time_bridge()
    print(f"  SYK saturates λ_L=2πT (MSS) ⟹ t* = (ln S)/(2πT); cosmic horizon T=H/2π ⟹")
    print(f"  t* = (ln S)/H = {br['t_star_over_inv_H']:.1f}/H  (Route-2 ln S≈{br['route2_lnS']:.0f})")
    out["diagnostics"]["scrambling_bridge"] = br

    # ── verdict ─────────────────────────────────────────────────────────────────
    rmt_ok = all(v["ok"] for v in gap.values())
    page_ok = page[12]["mean_ratio_q4"] > 0.97 and page[12]["mean_ratio_q2"] < page[12]["mean_ratio_q4"]
    otoc_ok = plateau > 1.8 and out["diagnostics"]["otoc"]["C_q4"][0] < 0.1
    out["verdict"] = dict(rmt_ok=rmt_ok, page_ok=page_ok, otoc_ok=otoc_ok,
                          all_ok=bool(rmt_ok and page_ok and otoc_ok))
    print("\n" + "=" * 74)
    print(f"VERDICT: SYK is maximally chaotic (RMT {rmt_ok}), its eigenstates are volume-law")
    print(f"(Page {page_ok}), and it scrambles area→volume completely (OTOC {otoc_ok}).")
    print("⟹ a maximally-scrambled (black-hole-dual) horizon HAS volume-law (Δ=1) entanglement —")
    print("  the physics SEDE's postulate invokes, now in the respected SYK model, not a toy.")
    print("  STILL OPEN (unchanged): that the COSMIC horizon is held in this driven steady")
    print("  state (Route 2 / the t*=(ln S)/H drive), and the magnitude. SYK gives the state,")
    print("  not the cosmology.")
    print("=" * 74)

    os.makedirs("results", exist_ok=True)
    with open("results/syk_scrambling.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote results/syk_scrambling.json")

    _make_figure(gap, page, ts, C4, C2)
    return out


def _make_figure(gap, page, ts, C4, C2):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(figure skipped: {e})")
        return
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    # (A) gap ratio
    Ns = sorted(gap)
    ax[0].plot(Ns, [gap[N]["r4"] for N in Ns], "o-", color="C3", label="SYK q=4 (chaotic)")
    ax[0].errorbar(Ns, [gap[N]["r4"] for N in Ns], yerr=[gap[N]["e4"] for N in Ns],
                   fmt="none", ecolor="C3", capsize=3)
    ax[0].plot(Ns, [gap[N]["r2"] for N in Ns], "s--", color="C0", label="q=2 (free/integrable)")
    for k, c in [("Poisson", "gray"), ("GOE", "C2"), ("GUE", "C1"), ("GSE", "C4")]:
        ax[0].axhline(syk.R_REF[k], ls=":", color=c, lw=1)
        ax[0].text(Ns[-1] + 0.18, syk.R_REF[k], k, va="center", ha="left", fontsize=8, color=c)
    ax[0].set_xticks(Ns)
    ax[0].set_xlim(Ns[0] - 0.4, Ns[-1] + 1.7)          # right margin so the RMT-class labels fit
    ax[0].set_xlabel("N (Majoranas)"); ax[0].set_ylabel(r"gap ratio $\langle r\rangle$")
    ax[0].set_title("(1) Maximal chaos: N mod 8 RMT class")
    ax[0].legend(loc="upper left", fontsize=8, framealpha=0.9)   # off the data, over the guide lines

    # (B) Page curve (N=12)
    N = 12
    p = page[N]
    ax[1].plot(p["n_A"], p["S_page"], "k--", label="Page (Haar, volume-law)")
    ax[1].plot(p["n_A"], p["S_q4"], "o-", color="C3", label="SYK q=4 eigenstates")
    ax[1].plot(p["n_A"], p["S_q2"], "s-", color="C0", label="q=2 free (sub-maximal)")
    ax[1].set_xlabel(r"subsystem size $n_A$ (qubits)"); ax[1].set_ylabel("entanglement S (nats)")
    ax[1].set_title("(2) Volume law: Page saturation (Δ=1)")
    ax[1].legend(fontsize=8)

    # (C) OTOC
    ax[2].plot(ts, C4, "-", color="C3", label="SYK q=4: full scramble (C→2)")
    ax[2].plot(ts, C2, "-", color="C0", label="q=2 free: incomplete")
    ax[2].axhline(2.0, ls=":", color="gray", lw=1)
    ax[2].set_xlabel("time t (1/J)"); ax[2].set_ylabel(r"OTOC $C(t)=2(1-\mathrm{Re}\,F)$")
    ax[2].set_title("(3) Dynamical area→volume scrambling")
    ax[2].legend(fontsize=8)

    fig.suptitle("SYK as the horizon scrambler — maximal chaos ⟹ volume-law (Δ=1) entanglement",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs("output", exist_ok=True)
    fig.savefig("results/fig_syk_scrambling.png", dpi=130)
    plt.close(fig)
    print("wrote results/fig_syk_scrambling.png")


if __name__ == "__main__":
    main()
