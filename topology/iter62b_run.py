"""
ITERATION 62b -- RECORDED REVISION of iter62_run.py after its method check FAILED: simulations made WITH an L = 20 Gpc cube were
not recognised (median Delta chi2 +12, i.e. the infinite model fitted them better), because 48 random cube orientations are too
coarse (typical misalignment ~25 deg). Fix: compute each torus covariance once in harmonic space and rotate it exactly with
real Wigner rotation matrices (C(R) = D C0 D^T), which allows 3000 orientations per L. Everything else unchanged. Larger
calibration: 50 isotropic sims, and cube sims at L = 20, 30, 35 Gpc (10 each) to measure what the test can and cannot see.
"""
import numpy as np, healpy as hp, time
from scipy.spatial.transform import Rotation
from scipy.special import eval_legendre
from iter62_lib import *
rng = np.random.default_rng(62)
F = "/mnt/user-data/uploads/COM_CMB_IQU-smica_1024_R2_02_full.fits"
I = hp.read_map(F, field=0)*1e6; TM = hp.read_map(F, field=3)
bl = hp.gauss_beam(np.radians(FWHM_DEG), lmax=LMAX_ISO); pl = np.ones(LMAX_ISO + 1)
alm = hp.map2alm(I, lmax=LMAX_ISO, iter=3)
dmap = hp.alm2map(hp.almxfl(alm, bl*pl), NSIDE, lmax=LMAX_ISO)
mfrac = hp.ud_grade(TM.astype(float), NSIDE); keep = np.where(mfrac > 0.9)[0]
th, ph = hp.pix2ang(NSIDE, keep); vec = np.array(hp.pix2vec(NSIDE, keep)).T; d = dmap[keep]; n = len(keep)
T = np.column_stack([np.ones(n), vec]); Qf, _ = np.linalg.qr(np.column_stack([T, rng.normal(size=(n, n - 4))])); Q = Qf[:, 4:]
dq = Q.T @ d
Yp = realY(th, ph, lm)*np.array([bl[l]*pl[l] for l, m in lm])[None, :]; YQ = Q.T @ Yp
cosg = np.clip(vec @ vec.T, -1, 1)
def iso_pix(lmin, lmax, cl):
    S = np.zeros((n, n))
    for l in range(lmin, lmax + 1): S += (2*l + 1)/(4*np.pi)*cl[l]*(bl[l]*pl[l])**2*eval_legendre(l, cosg)
    return S
SH_q = Q.T @ iso_pix(LMAX_T + 1, LMAX_ISO, cl_camb) @ Q
Sinf_q = Q.T @ iso_pix(2, LMAX_ISO, cl_camb) @ Q
KC = 35/CHI; clt = iso_cl_from_transfer(KC); rem = np.diag([cl_camb[l] - clt[l] for l, m in lm])
eps = 1.0; Agrid = np.exp(np.linspace(np.log(0.2), np.log(5), 121))
# real Wigner matrices: Y(R x) = D(R) Y(x), by least squares on sample points
Xs = rng.normal(size=(600, 3)); Xs /= np.linalg.norm(Xs, axis=1)[:, None]
def Yof(X): return realY(np.arccos(np.clip(X[:, 2], -1, 1)), np.arctan2(X[:, 1], X[:, 0]), lm)
Y0 = Yof(Xs); Y0pinv = np.linalg.pinv(Y0)
def Dmat(R): return (Y0pinv @ Yof(R.apply(Xs))).T
def chi2min(Sq, dvecs):
    w, U = np.linalg.eigh(Sq); P = (U.T @ dvecs)**2
    best = np.full(dvecs.shape[1], np.inf)
    for A in Agrid:
        lam = A*w + eps; best = np.minimum(best, (P/lam[:, None]).sum(0) + np.log(lam).sum())
    return best
def draw(Sq, A, k):
    w, U = np.linalg.eigh(Sq); return U @ (np.sqrt(np.clip(A*w + eps, 0, None))[:, None]*rng.normal(size=(len(w), k)))
c_inf_data = chi2min(Sinf_q, dq[:, None])
# amplitude for sims: best infinite-space A for the data
w0, U0 = np.linalg.eigh(Sinf_q); P0 = (U0.T @ dq)**2
A0 = Agrid[np.argmin([(P0/(A*w0 + eps)).sum() + np.log(A*w0 + eps).sum() for A in Agrid])]
C0 = {}
def C0_of(L):
    if L not in C0: C0[L] = torus_cov(L, Rotation.identity(), KC, 1.0) + rem
    return C0[L]
def torus_Sq(L, R):
    D = Dmat(R); return YQ @ (D @ C0_of(L) @ D.T) @ YQ.T + SH_q
groups = [("data", dq[:, None]), ("null", draw(Sinf_q, A0, 50))]
for Lt in (20000., 30000., 35000.):
    groups.append((f"cube{int(Lt/1000)}", np.column_stack([draw(torus_Sq(Lt, Rotation.random(random_state=rng)), A0, 1) for _ in range(10)])))
Dall = np.column_stack([g[1] for g in groups]); idx = np.cumsum([0] + [g[1].shape[1] for g in groups])
c_inf = chi2min(Sinf_q, Dall)
Ls = [20000., 25000., 27740., 30000., 35000., 40000., 50000., 60000.]
NOR = 3000; ROT = Rotation.random(NOR, random_state=np.random.default_rng(6201))
res = {}; t0 = time.time()
for L in Ls:
    b = np.full(Dall.shape[1], np.inf)
    for i in range(NOR): b = np.minimum(b, chi2min(torus_Sq(L, ROT[i]), Dall))
    res[L] = b - c_inf; print(f"L = {L/1000:.2f} Gpc done ({time.time() - t0:.0f} s)", flush=True)
np.save("iter62b_dchi.npy", np.array([res[L] for L in Ls]))
out = ["ITERATION 62b: wrap-around cube vs Planck SMICA (revision: exact rotations, 3000 orientations per L)", "",
       f"pixels {n} (Nside {NSIDE}, {FWHM_DEG:.0f} deg, TMASK > 0.9); amplitude for sims {A0:.2f}; Delta chi2 = torus(best orientation) - infinite",
       "", " L [Gpc]   DATA     null: median  5th pct   | cube sims (median): L=20    L=30    L=35"]
for L in Ls:
    r = res[L]; g = {nm: r[idx[i]:idx[i+1]] for i, (nm, _) in enumerate(groups)}
    out.append(f" {L/1000:6.2f}  {g['data'][0]:+7.2f}    {np.median(g['null']):+7.2f}  {np.percentile(g['null'], 5):+7.2f}    |"
               f"            {np.median(g['cube20']):+7.2f} {np.median(g['cube30']):+7.2f} {np.median(g['cube35']):+7.2f}")
# fraction of null sims at least as torus-like as the data, per L, and look-elsewhere over L
nulls = np.array([res[L][idx[1]:idx[2]] for L in Ls]); datav = np.array([res[L][0] for L in Ls])
out.append("")
out.append("chance of a null sky looking at least as cube-like as the data: " + "  ".join(f"{L/1000:.0f}: {np.mean(nulls[i] <= datav[i]):.2f}" for i, L in enumerate(Ls)))
best_null = (nulls - np.median(nulls, axis=1, keepdims=True)).min(0); best_data = (datav - np.median(nulls, axis=1)).min()
out.append(f"look-elsewhere over all L (data's most cube-like L relative to null median): p = {np.mean(best_null <= best_data):.2f}")
txt = "\n".join(out); print(txt); open("iter62b_result.txt", "w").write(txt + "\n")
