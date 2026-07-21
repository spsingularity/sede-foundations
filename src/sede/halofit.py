"""Halofit (Takahashi et al. 2012) nonlinear matter power spectrum. Δ²_nl from Δ²_lin(k)."""
import numpy as np


def _sigma2_gauss(R, k, D2lin):
    """σ²(R) with a Gaussian filter: ∫ dlnk Δ²_lin exp(-(kR)²)."""
    return np.trapezoid(D2lin * np.exp(-(k * R)**2), np.log(k))


def nonlinear_D2(k, D2lin, Om_z):
    """Nonlinear dimensionless power Δ²_nl(k) given linear Δ²_lin(k)=k³P/(2π²) and Ω_m(z).
    Takahashi 2012 calibration (flat, ~ΛCDM coefficients; adequate for a Fisher forecast)."""
    lnk = np.log(k)
    # find k_sigma where σ²(R)=1 (Gaussian filter)
    Rgrid = np.logspace(-2.5, 1.5, 400)
    s2 = np.array([_sigma2_gauss(R, k, D2lin) for R in Rgrid])
    if s2.min() > 1 or s2.max() < 1:
        return D2lin                      # no nonlinear scale in range → return linear
    lnR = np.log(Rgrid); lns2 = np.log(s2)
    Rs = np.exp(np.interp(0.0, lns2[::-1], lnR[::-1]))      # σ²(Rs)=1
    ksig = 1.0 / Rs
    # n_eff and C from derivatives of lnσ² at Rs
    d1 = np.gradient(lns2, lnR); d2 = np.gradient(d1, lnR)
    neff = -3.0 - np.interp(np.log(Rs), lnR, d1)
    Curv = -np.interp(np.log(Rs), lnR, d2)
    n = neff; om = Om_z
    # Takahashi 2012 coefficients
    an = 10**(1.5222 + 2.8553*n + 2.3706*n**2 + 0.9903*n**3 + 0.2250*n**4 - 0.6038*Curv)
    bn = 10**(-0.5642 + 0.5864*n + 0.5716*n**2 - 1.5474*Curv)
    cn = 10**(0.3698 + 2.0404*n + 0.8161*n**2 + 0.5869*Curv)
    gamman = 0.1971 - 0.0843*n + 0.8460*Curv
    alphan = abs(6.0835 + 1.3373*n - 0.1959*n**2 - 5.5274*Curv)
    betan = 2.0379 - 0.7354*n + 0.3157*n**2 + 1.2490*n**3 + 0.3980*n**4 - 0.1682*Curv
    mun = 0.0
    nun = 10**(5.2105 + 3.6902*n)
    y = k / ksig
    f = y / 4.0 + y**2 / 8.0
    D2Q = D2lin * ((1 + D2lin)**betan / (1 + alphan * D2lin)) * np.exp(-f)
    D2Hp = an * y**(3*1.0) / (1 + bn * y + (cn * y)**(3 - gamman))   # f1≈1 for flat
    D2H = D2Hp / (1 + mun / y + nun / y**2)
    return D2Q + D2H
