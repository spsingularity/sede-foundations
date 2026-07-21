"""Diagnostic plots + perturbation summary for the running-gamma SEDE-H model."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .gamma_computation import build_sigma_S_interp, gamma_eff_running
from .friedmann import E_SEDE_H_run, E_LCDM, compute_growth_factor
from .perturbations import sedeH_perturbation_report


def make_plots(Om=0.3185, s8=0.679, out='output/sedeH_running.png'):
    Ss = build_sigma_S_interp(p=1.0)
    D = lambda zz: compute_growth_factor(np.asarray(zz, float), Om)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))

    # (1) running gamma_eff(z) for p=1 vs the old constant gamma=1.5
    zr = np.linspace(0, 2.5, 30)
    gr = gamma_eff_running(D, p=1.0, z_arr=zr)['gamma_eff']
    ax[0].plot(zr, gr, 'b-', lw=2, label='p=1 extensive (Thm 4B)')
    ax[0].axhline(1.4964, color='r', ls='--', label='old constant γ=1.50 (p=5/3, z=0)')
    ax[0].axhline(1.0, color='g', ls=':', label='γ_data (const-γ MCMC)')
    ax[0].set_xlabel('z'); ax[0].set_ylabel(r'$\gamma_{\rm eff}(z)$')
    ax[0].set_title('Running structural coupling'); ax[0].legend(fontsize=8)

    # (2) f_sat(z): running p=1 vs constant gamma
    zz = np.linspace(0, 5, 200)
    Dz = D(zz)
    fsat_run = np.array([float(Ss(s8 * d)) / float(Ss(s8)) for d in Dz])
    g0 = 1.0
    fsat_const = (1 - np.exp(-g0 * Dz**2)) / (1 - np.exp(-g0))
    ax[1].plot(zz, fsat_run, 'b-', lw=2, label='running p=1')
    ax[1].plot(zz, fsat_const, 'r--', lw=1.5, label='constant γ=1.0')
    ax[1].set_xlabel('z'); ax[1].set_ylabel(r'$f_{\rm sat}(z)$')
    ax[1].set_title('Saturation function'); ax[1].legend(fontsize=8)

    # (3) EOS w(z): running SEDE-H vs DESI CPL band
    zw = np.linspace(0.001, 2.0, 250)
    Ef = E_SEDE_H_run(zw, Om, s8, Ss)
    rho = Ef**2 - Om * (1 + zw)**3; rho /= rho[0]
    w = -1 + (1/3.) * (1 + zw) * np.gradient(np.log(np.abs(rho)), zw)
    ax[2].plot(zw, w, 'b-', lw=2, label='running SEDE-H')
    w_desi = -0.838 - 0.62 * zw / (1 + zw)
    ax[2].plot(zw, w_desi, 'k--', lw=1.2, label='DESI CPL (w0=-0.84,wa=-0.62)')
    ax[2].axhline(-1, color='gray', ls=':')
    ax[2].set_xlabel('z'); ax[2].set_ylabel('w(z)'); ax[2].set_ylim(-1.3, -0.5)
    ax[2].set_title('Dark-energy EOS'); ax[2].legend(fontsize=8)

    fig.tight_layout(); fig.savefig(out, dpi=110)
    print('saved', out)

    print('\nRunning-SEDE-H best-fit perturbation observables:')
    sedeH_perturbation_report(Om, gamma=1.0, sigma8=s8,
                              E_run=lambda zz: E_SEDE_H_run(zz, Om, s8, Ss),
                              label='SEDE-H-running')


if __name__ == '__main__':
    make_plots()
