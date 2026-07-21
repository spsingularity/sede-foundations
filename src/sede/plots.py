"""
SEDE diagnostic plots.
Produces: E(z), f_sat(z), w_DE(z), posterior contours, cosmic history.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

from .friedmann import E_SEDE, E_LCDM, compute_growth_factor
from .theory import w_DE_algebraic, fsat_late

plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 13, 'legend.fontsize': 11,
    'figure.dpi': 150, 'lines.linewidth': 2,
})
COLORS = {'sede': '#D62728', 'lcdm': '#1F77B4', 'desi': '#2CA02C'}

def plot_E_fsat(Omega_m=0.302, gamma=1.4964, outdir='output/plots'):
    os.makedirs(outdir, exist_ok=True)
    z = np.linspace(0, 3, 500)
    E_s = E_SEDE(z, Omega_m, gamma)
    E_l = E_LCDM(z, Omega_m)
    D   = compute_growth_factor(z, Omega_m)
    fsat = fsat_late(z, 0.811 * D, 0.811, gamma)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel 1: E(z)
    ax = axes[0]
    ax.plot(z, E_s, color=COLORS['sede'], label='SEDE')
    ax.plot(z, E_l, color=COLORS['lcdm'], ls='--', label='ΛCDM')
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('E(z) = H(z)/H₀')
    ax.set_title('Hubble parameter')
    ax.legend()
    ax.set_xlim(0, 3)

    # Panel 2: E_SEDE / E_LCDM
    ax = axes[1]
    ax.plot(z, E_s / E_l, color='purple')
    ax.axhline(1, ls=':', color='gray')
    ax.axhline(0.95, ls=':', color='gray', alpha=0.5)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('E_SEDE / E_ΛCDM')
    ax.set_title('SEDE/ΛCDM ratio')
    ax.set_xlim(0, 3)
    ax.set_ylim(0.93, 1.02)

    # Panel 3: f_sat(z)
    ax = axes[2]
    ax.plot(z, fsat, color=COLORS['sede'], label='f_sat(z)')
    ax.axhline(1.0, ls=':', color='gray')
    ax.fill_between(z, fsat, 1, alpha=0.2, color=COLORS['sede'])
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('f_sat(z) = Ω_DE(z) / Ω_DE(0)')
    ax.set_title('Horizon saturation function')
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1.05)
    ax.legend()

    fig.suptitle(f'SEDE: Ω_m={Omega_m:.3f}, γ={gamma:.3f}', y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, 'sede_E_fsat.png')
    plt.savefig(path, bbox_inches='tight')
    plt.close()
    print(f'Saved: {path}')


def plot_w_eff(Omega_m_range=(0.28, 0.35), outdir='output/plots'):
    """Plot w_DE(0) as a function of Omega_m with DESI DR2 band."""
    os.makedirs(outdir, exist_ok=True)
    Om_arr = np.linspace(*Omega_m_range, 100)
    w_arr  = w_DE_algebraic(Om_arr)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(Om_arr, w_arr, color=COLORS['sede'], lw=2.5,
            label='SEDE: $w_{DE}(0) = (4\\Omega_m/3 - 1)/(1 - \\Omega_m)$')
    ax.axhline(-1, ls='--', color='gray', lw=1.5, label='Cosmological constant (Λ)')

    # DESI DR2 constraint: w_0 = -0.838 ± 0.055
    w_desi, sw_desi = -0.838, 0.055
    ax.axhspan(w_desi - sw_desi, w_desi + sw_desi,
               alpha=0.25, color=COLORS['desi'], label='DESI DR2 $w_0 = -0.838 \\pm 0.055$')
    ax.axhline(w_desi, color=COLORS['desi'], ls=':', lw=1)

    # Best-fit point
    Om_bf = 0.302
    w_bf  = w_DE_algebraic(Om_bf)
    ax.plot(Om_bf, w_bf, 'o', color=COLORS['sede'], ms=10,
            label=f'SEDE best fit: $\\Omega_m={Om_bf:.3f}$, $w={w_bf:.3f}$')

    ax.set_xlabel('$\\Omega_m$')
    ax.set_ylabel('$w_{DE}(0)$')
    ax.set_title('SEDE algebraic prediction vs DESI DR2 (Theorem 5)')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim(*Omega_m_range)
    ax.set_ylim(-1.05, -0.75)
    ax.grid(alpha=0.3)

    path = os.path.join(outdir, 'sede_w_algebraic.png')
    plt.savefig(path, bbox_inches='tight')
    plt.close()
    print(f'Saved: {path}')


def plot_cosmic_history(gamma=1.4964, Omega_m=0.302, outdir='output/plots'):
    """Plot Omega_DE(z) = f_sat(z) * Omega_DE0 across cosmic time."""
    os.makedirs(outdir, exist_ok=True)
    z = np.logspace(-2, 3.5, 1000)
    D   = compute_growth_factor(z, Omega_m)
    fsat = fsat_late(z, 0.811 * D, 0.811, gamma)
    OmDE0 = 1 - Omega_m - 9e-5
    OmDE  = fsat * OmDE0

    # Reheating: fsat goes from 1 to 0 (schematic)
    z_rh = np.linspace(1e10, 1e14, 100)
    Gamma_phi = 1e-2  # schematic decay rate
    fsat_rh = np.exp(-Gamma_phi / z_rh)
    fsat_rh /= fsat_rh.max()

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.semilogx(1 + z, OmDE, color=COLORS['sede'], lw=2, label='$\\Omega_{DE}(z)$ — SEDE')
    ax.axhline(OmDE0, ls=':', color='gray', lw=1.5, label=f'$\\Omega_{{DE,0}} = {OmDE0:.3f}$')
    ax.axhline(0.0, ls=':', color='gray', lw=1)
    ax.axvline(1 + 1089.9, ls='--', color='orange', lw=1.5, alpha=0.7, label='Recombination z=1090')
    ax.axvline(1 + 0.0,    ls='--', color='green', lw=1.5, alpha=0.7, label='Today z=0')

    ax.fill_between(1 + z, OmDE, alpha=0.2, color=COLORS['sede'])
    ax.set_xlabel('$1 + z$')
    ax.set_ylabel('$\\Omega_{DE}(z) = f_{sat}(z) \\cdot \\Omega_{DE,0}$')
    ax.set_title('SEDE Dark Energy History: $f_{sat}(z)$ traces structure growth')
    ax.legend(loc='upper right')
    ax.set_xlim(1, 1e4)
    ax.set_ylim(-0.02, 0.75)
    ax.grid(alpha=0.3)

    # Annotation
    ax.annotate('Matter era:\n$f_{sat} \\approx 0$', xy=(1e2, 0.05), fontsize=10,
                color='gray')
    ax.annotate('Structure era:\n$f_{sat}$ rising', xy=(3, 0.45), fontsize=10,
                color=COLORS['sede'])

    path = os.path.join(outdir, 'sede_cosmic_history.png')
    plt.savefig(path, bbox_inches='tight')
    plt.close()
    print(f'Saved: {path}')


def plot_posterior(chains_path='output/sede_chains.npy', outdir='output/plots'):
    """Plot 1D and 2D posteriors from MCMC chains."""
    os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(chains_path):
        print(f'[WARN] Chains not found at {chains_path}')
        return
    flat = np.load(chains_path)
    pnames = ['$\\Omega_m$', '$H_0$ [km/s/Mpc]', '$r_d$ [Mpc]']

    try:
        import corner
        fig = corner.corner(flat, labels=pnames,
                            quantiles=[0.16, 0.50, 0.84],
                            show_titles=True, title_fmt='.4f',
                            smooth=1.5, color=COLORS['sede'])
        path = os.path.join(outdir, 'sede_corner.png')
        fig.savefig(path, bbox_inches='tight')
        plt.close()
        print(f'Saved: {path}')
    except Exception as e:
        print(f'[WARN] Corner plot failed: {e}')


def plot_reheating_master_eq(outdir='output/plots'):
    """Plot f_sat evolution during reheating (Theorem 2)."""
    os.makedirs(outdir, exist_ok=True)
    from .theory import reheating_fsat

    # Schematic: H(t) = H_rh * (t_rh/t)^{1/2} for radiation-dominated
    # At reheating: H_rh = Gamma_phi
    Gamma_phi = 1.0  # units t_rh = 1
    H_rh = Gamma_phi

    def H_func(t):
        return H_rh / np.sqrt(t) if t > 0 else H_rh

    t_span = (0.01, 20.0)
    t, fsat = reheating_fsat(t_span, H_func, Gamma_phi, f0=1.0)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(t, fsat, color=COLORS['sede'], lw=2.5)
    ax.axhline(0, ls=':', color='gray')
    ax.axhline(1, ls=':', color='gray')
    ax.axvline(1.0, ls='--', color='orange', label='$t = \\Gamma_\\phi^{-1}$ (reheating)')
    ax.set_xlabel('$t / t_{rh}$')
    ax.set_ylabel('$f_{sat}(t)$ = inflaton fraction')
    ax.set_title('Reheating Master Equation (Theorem 2)\n'
                 '$\\dot{f}_{sat} = H f_{sat}(1-f_{sat}) - \\Gamma_\\phi f_{sat}$')
    ax.legend()
    ax.set_xlim(0, 15)
    ax.set_ylim(-0.05, 1.1)
    ax.grid(alpha=0.3)

    path = os.path.join(outdir, 'sede_reheating_fsat.png')
    plt.savefig(path, bbox_inches='tight')
    plt.close()
    print(f'Saved: {path}')


def make_all_plots():
    os.makedirs('output/plots', exist_ok=True)
    print('Generating SEDE diagnostic plots...')
    plot_E_fsat()
    plot_w_eff()
    plot_cosmic_history()
    plot_reheating_master_eq()
    plot_posterior()
    print('Done.')


if __name__ == '__main__':
    make_all_plots()
