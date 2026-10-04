"""Shared machinery for iteration 62 (wrap-around cube vs Planck SMICA, large scales)."""
import numpy as np, healpy as hp, camb
from scipy.special import sph_harm_y
from scipy.spatial.transform import Rotation
LMAX_T, LMAX_ISO, NSIDE, FWHM_DEG = 12, 30, 8, 20.0
CHI = 13870.0                                   # Mpc, comoving distance to last scattering (CAMB)
p = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.12, tau=0.054, As=2.1e-9, ns=0.965, lmax=80)
res = camb.get_results(p)
td = res.get_cmb_transfer_data(tp="scalar")
q = td.q; Ls = list(td.L)
T0 = 2.7255e6                                   # muK
def PR(k): return 2.1e-9*(k/0.05)**(0.965 - 1)
Delta = {l: td.delta_p_l_k[0, Ls.index(l), :] for l in range(2, LMAX_T + 1)}
def Delta_at(l, k): return np.interp(k, q, Delta[l], left=0, right=0)
cl_camb = res.get_unlensed_scalar_cls(CMB_unit="muK", raw_cl=True)[:, 0]
lm = [(l, m) for l in range(2, LMAX_T + 1) for m in range(-l, l + 1)]
LL = np.array([l for l, m in lm])
PHASE = np.cos(np.pi*(LL[:, None] - LL[None, :])/2)
def realY(theta, phi, lmax_set):
    out = []
    for l, m in lmax_set:
        Y = sph_harm_y(l, abs(m), theta, phi)
        if m > 0: out.append(np.sqrt(2)*(-1)**m*Y.real)
        elif m < 0: out.append(np.sqrt(2)*(-1)**m*Y.imag)
        else: out.append(Y.real)
    return np.array(out).T
def iso_cl_from_transfer(kcut):
    """C_l (muK^2) from the transfer functions, integral up to kcut -- normalisation check against CAMB."""
    out = {}
    for l in range(2, LMAX_T + 1):
        m = q < kcut; out[l] = 4*np.pi*np.trapezoid(PR(q[m])*Delta[l][m]**2/q[m], q[m])*T0**2
    return out
def torus_modes(Lmpc, kcut):
    nmax = int(np.ceil(kcut*Lmpc/(2*np.pi)))
    g = np.arange(-nmax, nmax + 1); n = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    n = n[np.any(n != 0, axis=1)]; k = 2*np.pi*n/Lmpc; kk = np.linalg.norm(k, axis=1); s = kk < kcut
    return k[s], kk[s]
def torus_cov(Lmpc, R, kcut, norm):
    k, kk = torus_modes(Lmpc, kcut); kh = (R.apply(k))/kk[:, None]
    th = np.arccos(np.clip(kh[:, 2], -1, 1)); ph = np.arctan2(kh[:, 1], kh[:, 0])
    Y = realY(th, ph, lm)                                                   # modes x coeffs
    Pk = (2*np.pi**2/kk**3)*PR(kk)/Lmpc**3
    D = np.stack([Delta_at(l, kk) for l, m in lm], 1)                      # modes x coeffs
    M = np.sqrt(Pk)[:, None]*D*Y
    return (4*np.pi)**2*(M.T @ M)*PHASE*T0**2*norm
