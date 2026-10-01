"""
Same-fashion test of four grid geometries (geometries.py) under two link rules:
  static   : every link keeps its conductance (1/length)
  adaptive : links thicken with the flux they carry (fungus / slime-mould rule), thickness s = min(sqrt(|Q|/q0), 1),
             conductance = s/length, iterated to steady state (q0 = 1 for all geometries)
Mass M = flux 4 pi M injected at the centre node; potential p = 0 on the outer shell. Measured from the shell-averaged
potential: g(r) = -dp/dr and v^2 = g r. Scores (same for every case):
  inner slope  d ln g / d ln r at r = 2-4   (Newton: -2)
  outer slope  d ln g / d ln r at r = 8-13  (MOND, flat rotation: -1; Newton -2; constant pull 0)
  Tully-Fisher slope of ln v_out^4 vs ln M   (MOND: 1; Newton-like outer: 2)
  anisotropy   spread of p across a shell at r = 10, in units of the potential drop per unit radius (0 = perfectly round)
  stable       did the adaptive iteration settle
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, sys, json
from scipy.sparse.csgraph import connected_components
from geometries import BUILDERS, R
def prepare(P, E):
    n = len(P); A = sp.coo_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(n, n))
    _, lab = connected_components(A, directed=False); keep = lab == lab[0]
    newid = -np.ones(n, int); newid[keep] = np.arange(keep.sum()); P = P[keep]
    E = E[keep[E[:, 0]] & keep[E[:, 1]]]; E = newid[E]
    Lng = np.linalg.norm(P[E[:, 0]] - P[E[:, 1]], axis=1)
    m = len(E); B = sp.csr_matrix((np.r_[np.ones(m), -np.ones(m)], (np.r_[np.arange(m), np.arange(m)], np.r_[E[:, 0], E[:, 1]])),
                                  shape=(m, len(P)))
    r = np.linalg.norm(P, axis=1); bnd = r > R - 1.2; bnd[0] = False
    return P, E, Lng, B, r, bnd
def solve(B, cond, free, n, M, x0=None):
    Lap = (B.T @ sp.diags(cond) @ B).tocsr(); A = Lap[free][:, free]
    s = np.zeros(n); s[0] = 4*np.pi*M
    x, _ = spl.cg(A, s[free], x0=x0, rtol=1e-9, maxiter=40000, M=sp.diags(1/A.diagonal()))
    p = np.zeros(n); p[free] = x; return p, x
def profile(p, r):
    bins = np.arange(1.0, R - 1.5, 1.0); mid = 0.5*(bins[1:] + bins[:-1])
    pm = np.array([p[(r >= a) & (r < b)].mean() for a, b in zip(bins[:-1], bins[1:])])
    g = -np.gradient(pm, mid); return mid, g
def run(geom, rule, M, seed=0, it=60):
    P, E, Lng, B, r, bnd = prepare(*BUILDERS[geom](seed)); n = len(P); free = ~bnd
    s_thk = np.ones(len(E)); x = None; settled = True
    for k in range(it if rule == "adaptive" else 1):
        p, x = solve(B, s_thk/Lng, free, n, M, x)
        if rule == "static": break
        Q = np.abs((s_thk/Lng)*(B @ p)); new = np.clip(np.sqrt(Q), 1e-4, 1.0)
        ch = np.max(np.abs(new - s_thk)); s_thk = 0.5*s_thk + 0.5*new
        if ch < 1e-4: break
    else:
        settled = rule == "static" or ch < 1e-3
    mid, g = profile(p, r)
    sh = (r > 9.5) & (r < 10.5); gi = np.interp(10, mid, g)
    aniso = np.std(p[sh])/max(abs(gi), 1e-12)
    return mid, g, aniso, settled
def slope(mid, g, lo, hi):
    m = (mid >= lo) & (mid <= hi) & (g > 0)
    return np.polyfit(np.log(mid[m]), np.log(g[m]), 1)[0] if m.sum() >= 2 else np.nan
if __name__ == "__main__":
    Ms = (4, 8, 16, 32); out = {}
    for geom in ("cubic", "mosaic", "web", "fungal"):
        for rule in ("static", "adaptive"):
            for seed in ((0, 1) if geom in ("mosaic", "fungal") else (0,)):
                rows = []; v4 = []
                for M in Ms:
                    mid, g, an, st = run(geom, rule, M, seed)
                    rows.append((slope(mid, g, 2, 4), slope(mid, g, 8, 13), an, st))
                    om = (mid >= 8) & (mid <= 13); v4.append(np.median((g*mid)[om])**2)
                tf = np.polyfit(np.log(Ms), np.log(v4), 1)[0]
                rr = np.array([x[:3] for x in rows]); stab = all(x[3] for x in rows)
                key = f"{geom}/{rule}/seed{seed}"
                out[key] = dict(inner=rr[:, 0].tolist(), outer=rr[:, 1].tolist(), aniso=rr[:, 2].tolist(), tf=tf, settled=stab)
                print(f"{key:24s} inner slope {np.nanmean(rr[:,0]):5.2f} | outer slope {np.nanmean(rr[:,1]):5.2f} "
                      f"(per M: {' '.join(f'{v:5.2f}' for v in rr[:,1])}) | TF slope {tf:4.2f} | anisotropy {np.mean(rr[:,2]):5.2f} | settled {stab}")
                sys.stdout.flush()
    json.dump(out, open("compare_geometries.json", "w"), indent=1)
