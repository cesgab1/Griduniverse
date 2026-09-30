"""
Separate true confinement from ordinary Coulomb falloff on small loops.
Free (Gaussian) field on the SAME grid and SAME loops:  q_i = J_i . M^+ . J_i,  M = B^T B (B = triangle-link incidence).
  U(1):  <W>_free = exp(-q/(2 beta))        SU(2) (leading order):  <W>_free = exp(-3 q/(2 beta))
Fit the Monte Carlo loops (per size class):  <W>_MC = e^{c0} * mean_i exp(-kappa * <W>_free exponent_i) * exp(-sigma * Area)
  kappa ~ 1  : coupling renormalisation;   sigma > 0 : genuine area law = confinement;   sigma = 0 : Coulomb phase.
"""
import numpy as np, sys, os, json, glob
os.environ["BOXL"] = "8"; sys.argv = ["x", "SU2", "1.0", "1"]
src = open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/strong_on_grid.py").read().split("group, beta = sys.argv")[0]
exec(src)
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg
from scipy.optimize import least_squares
Bm = csr_matrix((TD.ravel().astype(float), (np.repeat(np.arange(nt), 3), TL.ravel())), (nt, nl))
M = (Bm.T @ Bm).tocsr()
Jall = np.zeros((nl, len(loops)))
DISP = np.array([l[2] for l in links]); wind = np.zeros(len(loops), bool)
for k, (area, per, li, di) in enumerate(loops):
    wind[k] = np.linalg.norm((di[:, None] * DISP[li]).sum(0)) > 1e-6        # loop wraps around the periodic box
    if not wind[k]: np.add.at(Jall[:, k], li, di.astype(float))
_ar = np.array([l[0] for l in loops]); print("loops wrapping the box, by area:", {float(a): int(wind[np.isclose(_ar, a)].sum()) for a in np.unique(_ar)}, flush=True)
def block_cg(Bk, tol=1e-6, maxit=3000):
    X = np.zeros_like(Bk); b2 = (Bk * Bk).sum(0); act = b2 > 0          # empty loops (path doubles back on itself) have q = 0
    R = Bk.copy(); Pp = R.copy(); rr = b2.copy()
    for it in range(maxit):
        a = np.where(act)[0]
        if len(a) == 0: break
        AP = M @ Pp[:, a]; al = rr[a] / (Pp[:, a] * AP).sum(0)
        X[:, a] += Pp[:, a] * al; R[:, a] -= AP * al; rn = (R[:, a] ** 2).sum(0)
        Pp[:, a] = R[:, a] + Pp[:, a] * (rn / rr[a]); rr[a] = rn
        act[a[rn < tol**2 * b2[a]]] = False
    return X, it
QC = "/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/q_cache.npy"
q = list(np.load(QC)) if os.path.exists(QC) else []
for s0 in (range(0, len(loops), 240) if not q else []):
    X, it = block_cg(Jall[:, s0:s0+240]); q += list((Jall[:, s0:s0+240] * X).sum(0)); print("block", s0, "iterations", it, flush=True)
q = np.array(q); np.save(QC, q); print(f"free-field loop sums q computed: range {q.min():.2f}-{q.max():.2f}", flush=True)
areas = np.array([l[0] for l in loops]); pers = np.array([l[1] for l in loops])
out = {}
for f in sorted(glob.glob("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/strong_*_*.json")):
    d = json.load(open(f)); g, beta = d["group"], d["beta"]; cfac = 0.5 if g == "U1" else 1.5
    cls_q = []; y = []; e = []; A = []
    for R, Tt, w, we, _ in d["loops"]:
        m = np.isclose(areas, R * Tt) & np.isclose(pers, 2 * (R + Tt))
        cls_q.append(q[m]); y.append(w); e.append(we); A.append(R * Tt)
    y, e, A = map(np.array, (y, e, A)); ok = (y > 3 * e) & np.array([not wind[np.isclose(areas, R*Tt) & np.isclose(pers, 2*(R+Tt))].any() for R, Tt, *_ in d['loops']])
    def model(p, qs=cls_q):
        c0, kap, sig = p
        return np.array([c0 + np.log(np.mean(np.exp(-kap * cfac * qq / beta))) - sig * a for qq, a in zip(qs, A)])
    res = lambda p: ((model(p) - np.log(np.clip(y, 1e-9, None))) * y / e)[ok]
    fit = least_squares(res, [0, 1, 0]); J = fit.jac; cov = np.linalg.inv(J.T @ J) * max(1, (fit.fun**2).sum() / max(1, ok.sum() - 3))
    err = np.sqrt(np.diag(cov)); chi2 = (fit.fun**2).sum()
    free_only = least_squares(lambda p: res([p[0], p[1], 0.0]), [0, 1])
    print(f"{g:4s} beta={beta:3.1f}  plaquette {d['plaq']:.3f} | area-law term sigma = {fit.x[2]:+.4f} ± {err[2]:.4f} ({fit.x[2]/err[2]:+.1f} sd) | "
          f"coupling factor kappa = {fit.x[1]:.2f} ± {err[1]:.2f} | chi2 {chi2:.1f} / {ok.sum()-3} dof (free field alone: {(free_only.fun**2).sum():.1f})", flush=True)
    out[f"{g}_{beta}"] = dict(sigma=fit.x[2], err=err[2], kappa=fit.x[1], chi2=chi2, dof=int(ok.sum() - 3), chi2_free=(free_only.fun**2).sum(), plaq=d["plaq"])
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/strong_analysis.json", "w"), indent=1)
