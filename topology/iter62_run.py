"""ITERATION 62 (pre-registered in PREREG_62.md): wrap-around cube vs Planck SMICA on the largest scales."""
import numpy as np, healpy as hp, time
from scipy.spatial.transform import Rotation
from scipy.special import eval_legendre
from iter62_lib import *
rng = np.random.default_rng(62)
F = "/mnt/user-data/uploads/COM_CMB_IQU-smica_1024_R2_02_full.fits"
I = hp.read_map(F, field=0) * 1e6                                          # muK, RING
TM = hp.read_map(F, field=3)
bl = hp.gauss_beam(np.radians(FWHM_DEG), lmax=LMAX_ISO); pl = np.ones(LMAX_ISO + 1)   # pixel window not downloadable here; data and covariance use the same filter, so consistent
alm = hp.map2alm(I, lmax=LMAX_ISO, iter=3)
dmap = hp.alm2map(hp.almxfl(alm, bl*pl), NSIDE, lmax=LMAX_ISO)
mfrac = hp.ud_grade(TM.astype(float), NSIDE); keep = np.where(mfrac > 0.9)[0]
th, ph = hp.pix2ang(NSIDE, keep); vec = np.array(hp.pix2vec(NSIDE, keep)).T
d = dmap[keep]; n = len(keep)
# monopole/dipole projection
T = np.column_stack([np.ones(n), vec]); Qf, _ = np.linalg.qr(np.column_stack([T, rng.normal(size=(n, n - 4))]))
Q = Qf[:, 4:]; dq = Q.T @ d
# harmonic -> pixel matrices
Yp = realY(th, ph, lm)*np.array([bl[l]*pl[l] for l, m in lm])[None, :]
cosg = np.clip(vec @ vec.T, -1, 1)
def iso_pix(lmin, lmax, cl):
    S = np.zeros((n, n))
    for l in range(lmin, lmax + 1): S += (2*l + 1)/(4*np.pi)*cl[l]*(bl[l]*pl[l])**2*eval_legendre(l, cosg)
    return S
S_high = iso_pix(LMAX_T + 1, LMAX_ISO, cl_camb)
S_inf = iso_pix(2, LMAX_ISO, cl_camb)
KC = 35/CHI; clt = iso_cl_from_transfer(KC)
rem = np.array([cl_camb[l] - clt[l] for l, m in lm])                      # isotropic remainder above the k cut
eps = 1.0
Agrid = np.exp(np.linspace(np.log(0.2), np.log(5), 161))
def prep(S):
    Sq = Q.T @ S @ Q; w, U = np.linalg.eigh(Sq); return w, U
def chi2min(w, U, dvecs):
    P = (U.T @ dvecs)**2                                                    # (n-4) x nvec
    best = np.full(dvecs.shape[1], np.inf); bestA = np.zeros(dvecs.shape[1])
    for A in Agrid:
        lam = A*w + eps; c = (P/lam[:, None]).sum(0) + np.log(lam).sum()
        sel = c < best; best[sel] = c[sel]; bestA[sel] = A
    return best, bestA
# infinite model fit to data
w0, U0 = prep(S_inf); c_inf, A_inf = chi2min(w0, U0, dq[:, None]); A0 = A_inf[0]
# simulations: 30 isotropic (at fitted amplitude), 10 torus L = 20 Gpc (fresh orientations)
def draw(S, A, k):
    w, U = prep(S); return U @ (np.sqrt(np.clip(A*w + eps, 0, None))[:, None]*rng.normal(size=(len(w), k)))
sims_null = draw(S_inf, A0, 30)
sims_t20 = np.column_stack([draw(Yp @ (torus_cov(20000., Rotation.random(random_state=rng), KC, 1.0) + np.diag(rem)) @ Yp.T + S_high, A0, 1) for _ in range(10)])
D = np.column_stack([dq, sims_null, sims_t20])
c_inf_all, _ = chi2min(w0, U0, D)
Ls = [20000., 25000., 27740., 30000., 35000., 40000., 50000., 60000.]
ROT = Rotation.random(48, random_state=np.random.default_rng(6200))
best = {L: np.full(D.shape[1], np.inf) for L in Ls}
t0 = time.time()
for L in Ls:
    for i in range(len(ROT)):
        S = Yp @ (torus_cov(L, ROT[i], KC, 1.0) + np.diag(rem)) @ Yp.T + S_high
        w, U = prep(S); c, _ = chi2min(w, U, D); best[L] = np.minimum(best[L], c)
    print(f"L = {L/1000:.2f} Gpc done ({time.time() - t0:.0f} s)", flush=True)
out = ["ITERATION 62: wrap-around cube vs Planck SMICA, largest scales (expectations pre-registered)", "",
       f"pixels used {n} of {12*NSIDE**2} (Nside {NSIDE}, {FWHM_DEG:.0f} deg smoothing, TMASK > 0.9); fitted amplitude (infinite) {A0:.2f}",
       "", " L [Gpc]   data Delta chi2 (torus - infinite, best of 48 orientations)   null sims: median / 95%   L=20 sims: median"]
for L in Ls:
    dchi = best[L] - c_inf_all
    out.append(f" {L/1000:6.2f}      {dchi[0]:+8.2f}                                         {np.median(dchi[1:31]):+7.2f} / {np.percentile(dchi[1:31], 5):+7.2f}"
               f"          {np.median(dchi[31:]):+8.2f}")
out.append("(negative = torus preferred; '95%' column is the 5th percentile of the null sims, i.e. how negative chance alone gets)")
txt = "\n".join(out); print(txt); open("iter62_result.txt", "w").write(txt + "\n")
np.save("iter62_dchi.npy", np.array([best[L] - c_inf_all for L in Ls]))
