"""
Weak point 4: the strong force on the random grid.
4D Euclidean random grid (Poisson points, density 1, periodic box L^4), links = Delaunay edges, plaquettes = Delaunay triangles.
Field: a matrix on every link.  SU(2) (colour force; SU(3) confines the same way) vs U(1) (electromagnetism, piece 16).
Action: S = beta * sum_triangles (1 - Re tr U_triangle / N).  Monte Carlo: heat bath (SU(2)) / Metropolis (U(1)),
links updated in colour classes that share no triangle.
Measurement: Wilson loops ~ rectangles R x T drawn as shortest grid paths.  <W> ~ exp(-sigma*Area - mu*Perimeter):
  sigma > 0  -> the energy between two static charges grows linearly with distance = confinement
  sigma = 0  -> perimeter law, the force falls off = Coulomb phase (massless photon)
"""
import numpy as np, sys, json, time
from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
import os
L, PAD = float(os.environ.get("BOXL", 6)), 2.2

def build():
    n = rng.poisson(L**4); P = rng.uniform(0, L, (n, 4))
    imgs, offs, homes = [P], [np.zeros((n, 4), int)], [np.arange(n)]
    for s in np.array(np.meshgrid(*[[-1, 0, 1]] * 4)).reshape(4, -1).T:
        if not s.any(): continue
        Q = P + s * L; m = np.all((Q > -PAD) & (Q < L + PAD), axis=1)
        imgs.append(Q[m]); offs.append(np.repeat(s[None], m.sum(), 0)); homes.append(np.arange(n)[m])
    X = np.concatenate(imgs); home = np.concatenate(homes)
    tri = Delaunay(X).simplices
    # triangles with centroid in the home box (each physical triangle exactly once)
    T = set(); E = {}
    for comb in [(0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 2, 3), (0, 2, 4), (0, 3, 4), (1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]:
        t = np.sort(tri[:, comb], axis=1); c = X[t].mean(1); m = np.all((c >= 0) & (c < L), axis=1)
        for row in map(tuple, np.unique(t[m], axis=0)): T.add(row)
    for comb in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        e = np.sort(tri[:, comb], axis=1); c = X[e].mean(1); m = np.all((c >= 0) & (c < L), axis=1)
        for a, b in np.unique(e[m], axis=0): E[(a, b)] = None
    # link identity: (home a, home b, offset of b relative to a)
    off = np.concatenate(offs); key2id = {}; links = []
    for a, b in E:
        k = (home[a], home[b], tuple(off[b] - off[a]))
        if k not in key2id:
            key2id[k] = len(links); key2id[(home[b], home[a], tuple(off[a] - off[b]))] = -len(links) - 1
            links.append((home[a], home[b], X[b] - X[a]))
    def lid(a, b):
        v = key2id[(home[a], home[b], tuple(off[b] - off[a]))]
        return (v, 1) if v >= 0 else (-v - 1, -1)
    tris = [[lid(a, b), lid(b, c), lid(c, a)] for a, b, c in T]
    return n, P, links, tris

t0 = time.time(); n, P, links, tris = build(); nl, nt = len(links), len(tris)
TL = np.array([[x[0] for x in t] for t in tris]); TD = np.array([[x[1] for x in t] for t in tris])
print(f"4D grid: {n} points, {nl} links ({2*nl/n:.1f} per point), {nt} triangles ({nt*3/nl:.1f} per link)  [{time.time()-t0:.0f}s]", flush=True)

# ---- colour classes: links sharing a triangle go in different classes ----
nbr = [set() for _ in range(nl)]
for t in TL:
    for i in t: nbr[i].update(t)
col = -np.ones(nl, int)
for i in np.argsort([-len(s) for s in nbr]):
    used = {col[j] for j in nbr[i]}; c = 0
    while c in used: c += 1
    col[i] = c
ncol = col.max() + 1
inc_t, inc_p = [[] for _ in range(ncol)], [[] for _ in range(ncol)]
for ti, t in enumerate(TL):
    for p in range(3): inc_t[col[t[p]]].append(ti); inc_p[col[t[p]]].append(p)
inc_t = [np.array(x) for x in inc_t]; inc_p = [np.array(x) for x in inc_p]
cls = [np.where(col == c)[0] for c in range(ncol)]
print(f"{ncol} update classes", flush=True)

# ---- SU(2) as unit quaternions ----
def qmul(a, b):
    return np.concatenate([(a[..., 0]*b[..., 0] - (a[..., 1:]*b[..., 1:]).sum(-1))[..., None],
                           a[..., 0:1]*b[..., 1:] + b[..., 0:1]*a[..., 1:] + np.cross(a[..., 1:], b[..., 1:])], -1)
def qconj(a): return a * np.array([1, -1, -1, -1])
def oriented(U, d): return np.where(d[..., None] > 0, U, qconj(U))

def staple_sum_su2(U, c):
    """for links in class c: sum over their triangles of R (the rest of the loop, arranged so tr(U_l R) = tr(loop))"""
    t, p = inc_t[c], inc_p[c]
    L0 = oriented(U[TL[t, (p + 1) % 3]], TD[t, (p + 1) % 3]); L1 = oriented(U[TL[t, (p + 2) % 3]], TD[t, (p + 2) % 3])
    R = qmul(L0, L1)
    d = TD[t, p]; R = np.where(d[:, None] > 0, R, qconj(R))            # link traversed backwards: tr(U^dag R) = tr(U R^dag)
    S = np.zeros((nl, 4)); np.add.at(S, TL[t, p], R); return S[cls[c]]

def heatbath_su2(U, beta):
    for c in range(ncol):
        V = staple_sum_su2(U, c); k = np.linalg.norm(V, axis=1); Vh = V / k[:, None]
        a = beta * k                                                    # weight exp(a x0) sqrt(1-x0^2)   (S = -beta/2 tr)
        x0 = np.empty(len(a)); todo = np.arange(len(a))
        while len(todo):
            u = rng.random(len(todo)); aa = a[todo]
            x = 1 + np.log(u + (1 - u) * np.exp(-2 * aa)) / aa
            ok = rng.random(len(todo)) < np.sqrt(np.clip(1 - x**2, 0, 1))
            x0[todo[ok]] = x[ok]; todo = todo[~ok]
        v = rng.normal(size=(len(a), 3)); v /= np.linalg.norm(v, axis=1)[:, None]
        X = np.concatenate([x0[:, None], np.sqrt(1 - x0**2)[:, None] * v], 1)
        U[cls[c]] = qmul(X, qconj(Vh))                                   # U V = X |V|  ->  U = X V^-1
    return U

def plaq_su2(U):
    Ls = [oriented(U[TL[:, i]], TD[:, i]) for i in range(3)]
    return qmul(qmul(Ls[0], Ls[1]), Ls[2])[:, 0].mean()

# ---- U(1) ----
def metro_u1(th, beta, hits=3, step=1.2):
    for c in range(ncol):
        t, p = inc_t[c], inc_p[c]
        rest = TD[t, (p+1) % 3]*th[TL[t, (p+1) % 3]] + TD[t, (p+2) % 3]*th[TL[t, (p+2) % 3]]
        own = TD[t, p]; idx = cls[c]; pos = np.searchsorted(idx, TL[t, p])
        for _ in range(hits):
            new = th[idx] + rng.uniform(-step, step, len(idx))
            dS = np.zeros(len(idx))
            np.add.at(dS, pos, -beta*(np.cos(own*new[pos] + rest) - np.cos(own*th[idx][pos] + rest)))
            acc = rng.random(len(idx)) < np.exp(-dS); th[idx[acc]] = new[acc]
    return th
def plaq_u1(th): return np.cos((TD * th[TL]).sum(1)).mean()

# ---- Wilson loops along shortest grid paths ----
W = np.array([np.linalg.norm(l[2]) for l in links]); A_, B_ = np.array([l[0] for l in links]), np.array([l[1] for l in links])
G = coo_matrix((np.r_[W, W], (np.r_[A_, B_], np.r_[B_, A_])), (n, n)).tocsr()
dist, pred = shortest_path(G, directed=False, return_predecessors=True)
pairs = {}
for i, (a, b, _) in enumerate(links):
    if (a, b) not in pairs or W[i] < W[abs(pairs[(a, b)][0])]: pairs[(a, b)] = (i, 1); pairs[(b, a)] = (i, -1)
def path(a, b):
    out = []
    while b != a:
        p = pred[a, b]; out.append(pairs[(p, b)]); b = p
    return out[::-1]
def nearest(x): return np.argmin(np.linalg.norm((P - x + L/2) % L - L/2, axis=1))
DISPL = np.array([l[2] for l in links])
loops = []
sizes = [(1.5, 1.5), (1.5, 2), (2, 2), (2, 2.5), (1.5, 3), (2.5, 2.5), (2, 3), (2.5, 3), (3, 3), (2, 4), (3, 3.5), (3.5, 3.5)]
for R, Tt in sizes:
    got = 0
    while got < 80:
        mu, nu = rng.choice(4, 2, replace=False); o = rng.uniform(0, L, 4)
        e1, e2 = np.eye(4)[mu] * R, np.eye(4)[nu] * Tt
        cs = [nearest(o % L), nearest((o + e1) % L), nearest((o + e1 + e2) % L), nearest((o + e2) % L)]
        pth = sum((path(cs[i], cs[(i+1) % 4]) for i in range(4)), [])
        li_, di_ = np.array([x[0] for x in pth]), np.array([x[1] for x in pth])
        if len(li_) == 0 or np.linalg.norm((di_[:, None] * DISPL[li_]).sum(0)) > 1e-6: continue   # skip loops that wrap the box
        loops.append((R * Tt, 2 * (R + Tt), li_, di_)); got += 1
print(f"{len(loops)} Wilson loops prepared  [{time.time()-t0:.0f}s]", flush=True)

def measure(Ufield, group):
    out = []
    for area, per, li, di in loops:
        if group == "U1": out.append(np.cos((di * Ufield[li]).sum()))
        else:
            M = np.array([1., 0, 0, 0])
            for q in oriented(Ufield[li], di): M = qmul(M, q)
            out.append(M[0])
    return np.array(out)

group, beta = sys.argv[1], float(sys.argv[2])
therm, meas, gap = [int(x) for x in (sys.argv[4].split(",") if len(sys.argv) > 4 else (60, 60, 2))]
if group == "SU2": U = np.tile([1., 0, 0, 0], (nl, 1)); step = lambda U: heatbath_su2(U, beta); plaq = plaq_su2
else: U = np.zeros(nl); step = lambda U: metro_u1(U, beta); plaq = plaq_u1
for s in range(therm): U = step(U)
print(f"thermalised: plaquette {plaq(U):.4f}  [{time.time()-t0:.0f}s]", flush=True)
data, pl = [], []
for s in range(meas):
    for _ in range(gap): U = step(U)
    data.append(measure(U, group)); pl.append(plaq(U))
data = np.array(data)
areas = np.array([l[0] for l in loops]); pers = np.array([l[1] for l in loops]); nlinks = np.array([len(l[2]) for l in loops])
res = []
for R, Tt in sizes:
    m = np.isclose(areas, R * Tt) & np.isclose(pers, 2 * (R + Tt))
    per_cfg = data[:, m].mean(1)
    res.append((R, Tt, float(per_cfg.mean()), float(per_cfg.std() / np.sqrt(len(per_cfg))), float(nlinks[m].mean())))
    print(f"  loop {R} x {Tt}: <W> = {per_cfg.mean():.5f} ± {per_cfg.std()/np.sqrt(len(per_cfg)):.5f}", flush=True)
# fit ln W = c - sigma*A - mu*P over loops with a clear signal
A = np.array([[1, -r[0]*r[1], -2*(r[0]+r[1])] for r in res]); y = np.array([r[2] for r in res]); e = np.array([r[3] for r in res])
ok = y > 3 * e
wts = (y[ok] / e[ok])
coef, *_ = np.linalg.lstsq(A[ok] * wts[:, None], np.log(y[ok]) * wts, rcond=None)
# bootstrap error on sigma
sig = []
for _ in range(300):
    yb = y + rng.normal(0, e); okb = ok & (yb > 0)
    cb, *_ = np.linalg.lstsq(A[okb] * (yb[okb]/e[okb])[:, None], np.log(yb[okb]) * (yb[okb]/e[okb]), rcond=None); sig.append(cb[1])
print(f"RESULT {group} beta={beta}: plaquette {np.mean(pl):.4f}; string tension sigma = {coef[1]:.4f} ± {np.std(sig):.4f} (per l_d^2); perimeter mu = {coef[2]:.4f}; loops used {ok.sum()}/{len(res)}", flush=True)
json.dump(dict(per_loop=data.mean(0).tolist(), per_loop_err=(data.std(0)/np.sqrt(len(data))).tolist(), group=group, beta=beta, plaq=float(np.mean(pl)), sigma=float(coef[1]), sigma_err=float(np.std(sig)), mu=float(coef[2]), loops=res),
          open(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/strong_{group}_{beta}.json", "w"), indent=1)
