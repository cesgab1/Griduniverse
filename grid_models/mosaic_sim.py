"""
Random-mosaic simulation of the drained grid (2D and 3D, periodic).
Links = edges of a random Delaunay mosaic (rigid random network), stiffness EA/L.
Drained: every link's rest length becomes L(1+eps0) (Poisson lengthening) -> all slack; links pull only (tension-only).
Load: stretch the patch by e (isotropic or along one axis), let nodes relax.
Refill rule: a drained link that goes taut (l >= L(1+eps0)) reseals -> rest length back to L. A sealed link that goes loose (l < L) drains again.
Iterate to self-consistency. Output: stress sigma(e) relative to the fully sealed network at the same e  ->  mu(e).
Nothing is chosen except eps0, which only sets the SCALE (checked by running two values).
"""
import numpy as np, sys, json
from scipy.spatial import Delaunay
from scipy.optimize import minimize
from itertools import product
rng = np.random.default_rng(1)

def build(N, d):
    u = rng.random((N, d))
    shifts = np.array(list(product([-1, 0, 1], repeat=d)))
    P = np.concatenate([u + s for s in shifts]); idx = np.tile(np.arange(N), len(shifts)); sh = np.repeat(shifts, N, axis=0)
    tri = Delaunay(P)
    E = set()
    for simp in tri.simplices:
        for a in range(d+1):
            for b in range(a+1, d+1):
                i, j = simp[a], simp[b]
                if np.all(sh[i] == 0) or np.all(sh[j] == 0):
                    ii, jj = idx[i], idx[j]; s = sh[j] - sh[i]
                    key = (ii, jj, tuple(s)) if (ii, jj) < (jj, ii) or ii != jj else None
                    if ii > jj: key = (jj, ii, tuple(-s))
                    E.add((min(ii,jj),) and key)
    E = [e for e in E if e is not None]
    I = np.array([e[0] for e in E]); J = np.array([e[1] for e in E]); S = np.array([e[2] for e in E], float)
    return u, I, J, S

def lengths(x, B, I, J, S):
    dv = (x[J] - x[I] + S) @ B.T
    return dv, np.linalg.norm(dv, axis=1)

def relax(u, B, I, J, S, R, d):
    N = len(u)
    def f(z):
        x = z.reshape(N, d); dv, l = lengths(x, B, I, J, S)
        st = np.maximum(l - R, 0)                       # tension only
        En = 0.5*np.sum(st**2/R)
        fl = (st/R/np.maximum(l, 1e-15))[:, None]*dv    # dE/d(dv)
        gx = fl @ B                                      # back to reduced coords
        g = np.zeros_like(x); np.add.at(g, J, gx); np.add.at(g, I, -gx)
        return En, g.ravel()
    r = minimize(f, u.ravel(), jac=True, method="L-BFGS-B", options=dict(maxiter=5000, gtol=1e-13, ftol=1e-16))
    return r.x.reshape(N, d)

def stress(u, B, I, J, S, R, d, mode):
    dv, l = lengths(u, B, I, J, S); t = np.maximum(l - R, 0)/R     # force (EA=1)
    n = dv/np.maximum(l, 1e-15)[:, None]
    sig = np.einsum("i,i,ia,ib->ab", t, l, n, n)
    return np.trace(sig)/d if mode == "iso" else sig[0, 0]

def run(d, N, eps0, mode, es):
    u0, I, J, S = build(N, d)
    L = lengths(u0, np.eye(d), I, J, S)[1]
    out = []; u = u0.copy(); sealed = np.zeros(len(L), bool)
    for e in es:
        B = np.eye(d)*(1+e) if mode == "iso" else np.diag([1+e] + [1]*(d-1))
        for it in range(200):
            R = np.where(sealed, L, L*(1+eps0))
            u = relax(u, B, I, J, S, R, d); l = lengths(u, B, I, J, S)[1]
            new = sealed.copy(); new[(~sealed) & (l >= L*(1+eps0)*(1+1e-9))] = True; new[sealed & (l < L*(1-1e-9))] = False
            if np.array_equal(new, sealed): break
            sealed = new
        sg = stress(u, B, I, J, S, np.where(sealed, L, L*(1+eps0)), d, mode)
        # fully sealed reference (linear network) at same e
        ur = relax(u0.copy(), B, I, J, S, L, d); sr = stress(ur, B, I, J, S, L, d, mode)
        out.append((e, sg/sr if sr > 0 else 0, sealed.mean()))
        print(f"  d={d} {mode} e/eps0={e/eps0:7.3f}  mu={out[-1][1]:.4f}  sealed={sealed.mean():.3f}", flush=True)
    return out

if __name__ == "__main__":
    d, N, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    res = {}
    for eps0 in (0.01, 0.02):
        es = eps0*np.concatenate([np.linspace(0.05, 0.9, 12), np.linspace(0.95, 1.2, 11), np.linspace(1.3, 4, 10)])
        res[eps0] = run(d, N, eps0, mode, es)
    json.dump({str(k): v for k, v in res.items()}, open(f"mosaic_{d}d_{mode}.json", "w"))
