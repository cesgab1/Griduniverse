"""
What if each grid point is a quantum tetrahedron (as in loop quantum gravity / group field theory)?
Tetrahedron vertices r_v = (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1): exactly the tri-bimaximal 'diagonal' directions.
Edge midpoints m_e = +-(1,0,0), +-(0,1,0), +-(0,0,1): exactly the 'axis' directions.
Fields live ON the building block:
   vertex field x_v (4 numbers)  -> its triplet part  phi = sum_v x_v r_v   (neutrino-sector flavon)
   edge field   y_e (6 numbers)  -> its triplet part  chi = sum_e y_e m_e   (charged-lepton flavon)
Interactions are LOCAL on the tetrahedron, all with random O(1) strengths:
   on-site      sum_v g(x_v) + sum_e h(y_e)                    (g, h quartic polynomials)
   incidence    sum over (vertex on edge) k(x_v, y_e)          (12 pairs)
   neighbours   sum over vertex pairs p(x_v, x_w) + adjacent edge pairs q(y_e, y_f)
The building block's symmetry makes these automatically tetrahedral-invariant; nothing about directions is put in by hand.
Measured: misalignment of phi from the nearest diagonal and of chi from the nearest axis at the global minimum.
Compared with the generic-potential results (A4 with all 21 invariants: 0% within 0.15 rad; cross-couplings suppressed: ~8%).
"""
import numpy as np, itertools, json, sys
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
R = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], float)
E = list(itertools.combinations(range(4), 2)); M = np.array([(R[a] + R[b]) / 2 for a, b in E])     # edge midpoints (axes)
inc = [(v, i) for i, (a, b) in enumerate(E) for v in (a, b)]
vpairs = list(itertools.combinations(range(4), 2))
epairs = [(i, j) for i, j in itertools.combinations(range(6), 2) if set(E[i]) & set(E[j])]        # edges sharing a vertex
DIAG = [r / np.sqrt(3) for r in R]; AXES = [np.eye(3)[i] for i in range(3)]
def ang(v, T):
    n = np.linalg.norm(v); return np.pi/2 if n < 1e-8 else min(np.arccos(min(1, abs(v @ t) / n)) for t in T)
def poly1(c, z): return c[0]*z + c[1]*z**2 + c[2]*z**3 + c[3]*z**4
def poly2(c, a, b):   # symmetric-enough generic 2-variable polynomial up to quartic (no pure terms; those are on-site)
    return c[0]*a*b + c[1]*(a*a*b) + c[2]*(a*b*b) + c[3]*a*a*b*b + c[4]*(a**3*b) + c[5]*(a*b**3)
def poly2s(c, a, b):  # symmetric under a<->b (same kind of object)
    return c[0]*a*b + c[1]*(a*a*b + a*b*b) + c[2]*a*a*b*b + c[3]*(a**3*b + a*b**3)
def make(scale_cross):
    g = rng.normal(0, 1, 4); g[1] = -abs(rng.normal(1, .3)); g[3] = abs(rng.normal(1, .3)) + 0.3
    h = rng.normal(0, 1, 4); h[1] = -abs(rng.normal(1, .3)); h[3] = abs(rng.normal(1, .3)) + 0.3
    k = rng.normal(0, 1, 6) * scale_cross; p = rng.normal(0, 1, 4) * 0.5; q = rng.normal(0, 1, 4) * 0.5
    def V(z):
        x, y = z[:4], z[4:]
        v = sum(poly1(g, xi) for xi in x) + sum(poly1(h, yi) for yi in y)
        v += sum(poly2(k, x[a], y[i]) for a, i in inc)
        v += sum(poly2s(p, x[a], x[b]) for a, b in vpairs) + sum(poly2s(q, y[i], y[j]) for i, j in epairs)
        return v
    return V
def bounded(V):
    d = rng.normal(size=(400, 10)); d /= np.linalg.norm(d, axis=1)[:, None]
    return all(V(20 * u) > V(np.zeros(10)) + 50 for u in d)
def run(scale_cross, n):
    mis, dead, kept = [], 0, 0
    while kept < n:
        V = make(scale_cross)
        if not bounded(V): continue
        kept += 1
        best = min((minimize(V, rng.normal(0, 1.5, 10), method="BFGS") for _ in range(20)), key=lambda r: r.fun)
        x, y = best.x[:4], best.x[4:]
        phi, chi = x @ R, y @ M
        if min(np.linalg.norm(phi), np.linalg.norm(chi)) < 0.05 * max(np.linalg.norm(phi), np.linalg.norm(chi), 1e-9): dead += 1; mis.append(np.pi/2); continue
        mis.append(max(ang(phi, DIAG), ang(chi, AXES)))
    mis = np.array(mis)
    return mis, dead
out = {}
for sc in [float(a) for a in (sys.argv[3].split(",") if len(sys.argv) > 3 else ["1.0", "0.3", "0.1"])]:
    mis, dead = run(sc, int(sys.argv[2]) if len(sys.argv) > 2 else 200)
    print(f"tetrahedral building block, vertex-edge coupling x{sc}: {len(mis)} random local potentials | a triplet field left unbroken in {dead/len(mis)*100:.0f}% | "
          f"aligned within 0.02 rad {np.mean(mis<0.02)*100:.1f}%, 0.15 rad {np.mean(mis<0.15)*100:.1f}%, 0.3 rad {np.mean(mis<0.3)*100:.1f}% | median {np.median(mis):.2f} rad", flush=True)
    out[sc] = dict(mis=mis.tolist(), dead=dead)
TAG = sys.argv[3] if len(sys.argv) > 3 else "all"
json.dump(out, open(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/tetra_{TAG}.json", "w"))
