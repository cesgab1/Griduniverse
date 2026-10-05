"""Iteration 100: fit each era separately (PREREG_100.md)."""
import os, numpy as np
from scipy.optimize import minimize_scalar
here = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(here, "iter99_timescape.py")).read().split("# self-check")[0])
full_icov, Z, M, ZH = sn_icov.copy(), z_sn.copy(), mb.copy(), zhel.copy()
def sn_slice(DM, m):
    ic = np.linalg.inv(np.linalg.inv(full_icov)[np.ix_(m, m)])
    r = M[m] - 5*np.log10((1 + ZH[m])*np.interp(Z[m], zgrid, DM)); B = (ic @ r).sum()
    return r @ ic @ r - B**2/ic.sum()
bao_all = bao.copy(); icov_all = np.linalg.inv(bao_icov)
def bao_slice(DM, m):
    global bao, bao_icov
    bao = bao_all[m].reset_index(drop=True); bao_icov = np.linalg.inv(icov_all[np.ix_(m, m)])
    c = bao_chi2(DM); bao, bao_icov = bao_all, np.linalg.inv(icov_all); return c
def fit(fun, lo, hi):
    r = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-5)); x0, c0 = r.x, r.fun
    grid = np.linspace(lo, hi, 400); cs = np.array([fun(g) for g in grid]); ok = grid[cs <= c0 + 1]
    return x0, (ok.min(), ok.max())
mk = {"timescape": (timescape, 0.30, 0.99), "LCDM": (lambda p: flat("LCDM", p), 0.05, 0.95)}
eras = []
if SNSET != "UNION3":
    for lo, hi in ((0.0, 0.2), (0.2, 0.5), (0.5, 9)):
        m = (Z >= lo) & (Z < hi); eras.append((f"SN z {lo}-{hi if hi < 9 else 'max'} (N={m.sum()})", "sn", m))
else:
    for lo, hi in ((0.0, 0.2), (0.2, 0.5), (0.5, 9)):
        m = (Z >= lo) & (Z < hi); eras.append((f"SN z {lo}-{hi if hi < 9 else 'max'} (N={m.sum()} bins)", "sn", m))
zb = bao_all.z.values
eras += [("BAO z < 1", "bao", zb < 1), ("BAO z > 1", "bao", zb > 1)]
out = [f"Iteration 100 -- SNSET={SNSET}   (best value [68% range])"]
for name, (f, lo, hi) in mk.items():
    out.append(f" {name}: {'f_v0' if name == 'timescape' else 'Om'} per era")
    vals = []
    for lab, kind, m in eras:
        fun = (lambda p, m=m: sn_slice(f(p), m)) if kind == "sn" else (lambda p, m=m: bao_slice(f(p), m))
        x, (a, b) = fit(fun, lo, hi); vals.append((x, a, b))
        edge = " (hits range edge)" if a <= lo + 1e-3 or b >= hi - 1e-3 else ""
        out.append(f"   {lab:28s} {x:.3f} [{a:.3f}-{b:.3f}]{edge}")
    # largest pairwise tension (symmetric half-widths)
    t = 0; pair = ""
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            xi, ai, bi = vals[i]; xj, aj, bj = vals[j]
            si, sj = (bi - ai)/2, (bj - aj)/2
            if si > 0 and sj > 0:
                tt = abs(xi - xj)/np.hypot(si, sj)
                if tt > t: t, pair = tt, f"{eras[i][0]} vs {eras[j][0]}"
    out.append(f"   largest disagreement between eras: {t:.1f} sigma ({pair})")
txt = "\n".join(out); print(txt); open(os.path.join(here, f"iter100_{SNSET}.txt"), "w").write(txt + "\n")
