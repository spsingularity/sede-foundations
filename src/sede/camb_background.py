"""CAMB-exact background for the Phase-1 joint (CAMB in the loop).

⚠️ DEPRECATED for headline use. The model='sedeH' background here is generated from
sede.friedmann.E_SEDE_H — the λ=1 "temperature-factor" COUSIN (w0≈−0.86, does NOT cross −1),
NOT the canonical Barrow λ=0.5 model (E_SEDE_lambda/E_SEDE_volume, w0≈−1.0, crosses −1). This
module is imported only by the deprecated run_camb_joint.py / run_camb_mcmc.py and is NOT in
reproduce_all.py. For the canonical headline joint use run_lambda_verify.py / run_barrow_mcmc.py
(ΔDIC≈−2.9) and run_joint_fullcmb.py (full-CMB joint, Δχ²=−3.17). See paper §3.2 for the
canonical-vs-cousin distinction.

Single engine for ALL geometry + sound horizons: SEDE-H enters CAMB as a tabulated
w(a) dark-energy fluid (PPF; c_s²=1, which Theorem 6 shows is the correct smooth-DE
treatment). CAMB then provides D_M(z), D_H(z), D_L(z), r_drag, r_star, z_star — so
BAO (with rd=r_drag, no free fudge), SN, and the compressed CMB (R, l_A) are all
mutually consistent and Planck-accurate (validated: LCDM@Planck → R −0.7σ, l_A +0.6σ).

w(a) is generated from sede.friedmann.E_SEDE_H (robust: it is a ratio, and radiation
is negligible at z<5 where DE matters); CAMB supplies the accurate radiation+neutrino
background for the absolute distances. Background-only calls are ~0.02–0.05 s.
"""
from __future__ import annotations

import numpy as np
import camb
from camb import dark_energy

from .friedmann import E_SEDE_H

C_KM_S = 299792.458


def _w_of_a_sedeH(Om, gamma, n=400, zmax=8.0):
    z = np.linspace(0.0, zmax, n)
    E = E_SEDE_H(z, Om, gamma)
    rho = np.maximum(E**2 - Om * (1 + z)**3, 1e-8)
    w = -1.0 + (1.0 / 3.0) * (1 + z) * np.gradient(np.log(rho), z)
    a = 1.0 / (1.0 + z)
    i = np.argsort(a)
    return a[i], w[i]


class Background:
    """CAMB background for model='lcdm' or 'sedeH' (SEDE-H via w(a) table)."""

    def __init__(self, Om, H0, gamma=1.0, ombh2=0.02237, model='sedeH', mnu=0.06):
        self.Om, self.H0, self.gamma, self.model = Om, H0, gamma, model
        h = H0 / 100.0
        omch2 = Om * h**2 - ombh2
        if omch2 <= 0:
            raise ValueError("omch2<=0")
        p = camb.CAMBparams()
        p.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=0.0)
        if model in ('sedeH', 'cousin'):   # naming-unification alias (physics identical)
            a, w = _w_of_a_sedeH(Om, gamma)
            de = dark_energy.DarkEnergyPPF()
            de.set_w_a_table(a, w)
            p.DarkEnergy = de
        p.WantCls = False
        p.WantTransfer = False
        self._r = camb.get_background(p)
        d = self._r.get_derived_params()
        self.r_drag = d['rdrag']
        self.r_star = d['rstar']
        self.z_star = d['zstar']
        self.theta_star = d['thetastar']        # 100*theta_star

    # distances (Mpc)
    def D_M(self, z):
        return self._r.comoving_radial_distance(np.atleast_1d(z))

    def D_H(self, z):
        return C_KM_S / self._r.hubble_parameter(np.atleast_1d(z))

    def D_V(self, z):
        z = np.atleast_1d(z); dm = self.D_M(z); dh = self.D_H(z)
        return (z * dm**2 * dh) ** (1.0 / 3.0)

    def D_L(self, z):
        z = np.atleast_1d(z)
        return (1.0 + z) * self.D_M(z)

    def distance_modulus(self, z):
        return 5.0 * np.log10(np.clip(self.D_L(z), 1e-8, None)) + 25.0

    def hubble(self, z):
        return self._r.hubble_parameter(np.atleast_1d(z))

    # compressed CMB
    def R_shift(self):
        DMz = float(self._r.comoving_radial_distance(self.z_star))
        return np.sqrt(self.Om) * self.H0 * DMz / C_KM_S

    def l_A(self):
        DMz = float(self._r.comoving_radial_distance(self.z_star))
        return np.pi * DMz / self.r_star
