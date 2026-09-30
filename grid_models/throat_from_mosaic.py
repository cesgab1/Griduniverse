"""
Throat strengths from the random mosaic itself (3D Poisson-Voronoi, periodic).
Fluid pockets = mosaic cells. A pocket keeps its fluid when the pull on its fluid (∝ g * volume V) beats the
void suction acting through its walls (∝ s * area A):   sealed  iff  g > s * (A/V - t_c)
   -> order in which pockets reseal = ascending A/V (big roomy cells hold fluid first). Shape of A/V spread: pure geometry.
Links = Delaunay edges between cell seeds; a link is sealed when BOTH its pockets are sealed (variant: EITHER).
Stiffness measured directly on the network (no linear assumption). Zero pull = marginal point (stiffness just 0).
Output: mu(x) = K(x)/K_full vs x ∝ (A/V - t_c); then SPARC (a0 scale only) and Cassini Q2.
"""
import numpy as np, json
from scipy.spatial import Voronoi, ConvexHull
from itertools import product
exec(open("drained_stiffness.py").read().split("res = {}")[0])      # build(), moduli()

def cells(u):
    N, d = u.shape; shifts = np.array(list(product([-1, 0, 1], repeat=3)))
    P = np.concatenate([u + s for s in shifts]); vor = Voronoi(P)
    c0 = int(np.where((shifts == 0).all(1))[0][0])*N
    AV = np.empty(N)
    for i in range(N):
        reg = vor.regions[vor.point_region[c0+i]]
        h = ConvexHull(vor.vertices[reg]); AV[i] = h.area/h.volume
    return AV

out = {}
for rule in ("both", "either"):
    curves = []
    for sd in range(3):
        rng = np.random.default_rng(200+sd); N = 1200
        u, I, J, n, L = build(N, 3, rng)
        AV = cells(u)*(1/N)**(1/3)                         # dimensionless (mean cell size = 1)
        cuts = np.quantile(AV, np.concatenate([np.linspace(0.05, 0.95, 46), 1-np.logspace(-1.35, -3, 8), [1.0]]))
        rows = []
        for t in cuts:
            s = AV <= t
            keep = (s[I] & s[J]) if rule == "both" else (s[I] | s[J])
            kb, ks = moduli(N, 3, I, J, n, L, keep)
            rows.append((t, keep.mean(), kb, ks))
        rows = np.array(rows); rows[:, 2] /= rows[-1, 2]; rows[:, 3] /= rows[-1, 3]
        curves.append(rows)
    C = np.mean(curves, 0); K = 0.5*(C[:, 2] + C[:, 3])
    i0 = max(np.argmax(K > 1e-3) - 1, 0); tc = C[i0, 0]                # marginal point: last cut with no stiffness
    x = C[:, 0] - tc
    out[rule] = dict(t=C[:, 0].tolist(), f=C[:, 1].tolist(), K=K.tolist(), tc=tc)
    print(f"\nrule '{rule}': marginal point at A/V = {tc:.3f} (sealed-link fraction {C[i0,1]:.3f})")
    print("   A/V - tc    sealed links   K/K_full")
    for j in list(range(i0, 46, 3)) + list(range(46, len(K))): print(f"   {x[j]:8.4f}     {C[j,1]:.3f}       {K[j]:.4f}")
    # deep-regime exponent: K ∝ x^p near the marginal point (MOND needs p = 1)
    m = (K > 0.02) & (K < 0.3) & (x > 0)
    p = np.polyfit(np.log(x[m]), np.log(K[m]), 1)[0]
    print(f"   near the marginal point K ∝ (pull)^{p:.2f}   (square-root law needs 1.00)")
    print(f"   tail: 1-mu at 2x and 3x the pull where mu=0.5: ", end="")
    xh = np.interp(0.5, K, x); print(f"{1-np.interp(2*xh, x, K):.4f}, {1-np.interp(3*xh, x, K):.4f}")
json.dump(out, open("throat_from_mosaic.json", "w"))
