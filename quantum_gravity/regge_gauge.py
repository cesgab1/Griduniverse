"""
Geometric version: gravity carried by the LINK LENGTHS of a random grid (Regge calculus, 1961; BORROWED).
3-D random Delaunay grid (the Mosaic's links). Action S = sum over links of length x deficit angle (curvature
lives on links as the angle missing around them). Question: is 'moving the nodes' an EXACT symmetry
(zero-cost directions of the action, 3 per interior node) on a random, irregular grid?
Hessian H_ef = d(deficit_e)/d(l_f) over interior links (Schlafli identity: dS/dl_e = deficit_e), by complex-step
differentiation (no subtraction error). Run with GRID=poisson (fully random) or GRID=jitter (randomly shaken lattice).
"""
import numpy as np, os
from scipy.spatial import Delaunay, ConvexHull
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(4)
MODE = os.environ.get('GRID','poisson')
if MODE=='poisson': X = rng.uniform(0, 1, (160, 3))
else:
    g=np.arange(6)/5.0; X=np.array(np.meshgrid(g,g,g,indexing='ij')).reshape(3,-1).T; X=X+rng.uniform(-0.35,0.35,X.shape)/5
tri = Delaunay(X); T = tri.simplices
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
E = np.unique(np.sort(np.vstack([T[:, list(p)] for p in pairs]), axis=1), axis=0); eid = {tuple(e): k for k, e in enumerate(E)}
hull = ConvexHull(X); bverts = set(np.unique(hull.simplices))
hull_e = set(eid[tuple(sorted((f[a], f[b])))] for f in hull.simplices for a, b in ((0, 1), (0, 2), (1, 2)))
interior = np.array([k for k in range(len(E)) if k not in hull_e]); ivert = np.array([v for v in range(len(X)) if v not in bverts])
TE = np.array([[eid[tuple(sorted((t[a], t[b])))] for a, b in pairs] for t in T])      # 6 edge ids per tet
def dihedrals(l):
    """dihedral angle of each tet at each of its 6 edges, from edge lengths only (intrinsic)."""
    L2 = l[TE]**2; th = np.zeros(L2.shape, dtype=L2.dtype)
    # vertices 0..3; edge index map for pair (i,j)
    pid = {p: k for k, p in enumerate(pairs)}; pid.update({(b, a): k for (a, b), k in pid.items()})
    for k, (a, b) in enumerate(pairs):
        c, d = [v for v in range(4) if v not in (a, b)]
        # vectors from a: u=b-a, v=c-a, w=d-a ; Gram from lengths
        uu = L2[:, pid[(a, b)]]; vv = L2[:, pid[(a, c)]]; ww = L2[:, pid[(a, d)]]
        uv = (uu + vv - L2[:, pid[(b, c)]])/2; uw = (uu + ww - L2[:, pid[(b, d)]])/2; vw = (vv + ww - L2[:, pid[(c, d)]])/2
        n1n2 = uu*vw - uw*uv; n1 = uu*vv - uv**2; n2 = uu*ww - uw**2
        th[:, k] = np.arccos(n1n2/np.sqrt(n1*n2)) if np.iscomplexobj(L2) else np.arccos(np.clip(n1n2/np.sqrt(n1*n2), -1, 1))
    return th
def deficits(l):
    s = np.zeros(len(E), dtype=l.dtype); np.add.at(s, TE.ravel(), dihedrals(l).ravel()); return 2*np.pi - s
l0 = np.linalg.norm(X[E[:, 0]] - X[E[:, 1]], axis=1); eps0 = deficits(l0)[interior]
H = np.zeros((len(interior), len(interior))); h = 1e-30
for c, f in enumerate(interior):
    q = l0.astype(complex); q[f] += 1j*h; H[:, c] = deficits(q)[interior].imag/h     # complex step: exact to rounding
asym = np.abs(H - H.T).max()/np.abs(H).max(); H = (H + H.T)/2; ev = np.linalg.eigvalsh(H)
# gauge directions: lengths change when interior nodes move, J = dl/dx
J = np.zeros((len(interior), 3*len(ivert))); pos = {e: r for r, e in enumerate(interior)}
for c, v in enumerate(ivert):
    for k in np.where((E[:, 0] == v) | (E[:, 1] == v))[0]:
        if k not in pos: continue
        o = E[k, 1] if E[k, 0] == v else E[k, 0]; J[pos[k], 3*c:3*c + 3] = (X[v] - X[o])/l0[k]
gres = np.linalg.norm(H @ J)/(np.linalg.norm(H)*np.linalg.norm(J))
rnd = rng.normal(size=J.shape); rres = np.linalg.norm(H @ rnd)/(np.linalg.norm(H)*np.linalg.norm(rnd))
scale = np.abs(ev).max(); nz = int(np.sum(np.abs(ev) < 1e-6*scale))
out = [f"random 3-D grid: {len(X)} nodes ({len(ivert)} interior), {len(E)} links ({len(interior)} interior), {len(T)} tetrahedra",
       f"flat start: max |deficit| on interior links = {np.abs(eps0).max():.1e} (flat space = 0)",
       f"Hessian symmetry check (should be ~0): {asym:.1e}",
       f"zero-cost directions found: {nz};  expected from moving interior nodes: 3 x {len(ivert)} = {3*len(ivert)}",
       f"gauge test |H.J|/(|H||J|) = {gres:.1e}   (random directions for comparison: {rres:.1e})",
       f"other eigenvalues: {np.sum(ev > 1e-6*scale)} positive, {np.sum(ev < -1e-6*scale)} negative (3-D gravity has no waves: all others are non-zero)"]
out.insert(0, f"grid: {MODE}; complex-step derivative (exact to rounding); min/max dihedral {np.degrees(dihedrals(l0).min()):.2f}/{np.degrees(dihedrals(l0).max()):.2f} deg"); txt = "\n".join(out); print(txt); open(os.path.join(HERE, f"regge_gauge_{MODE}.txt"), "w").write(txt + "\n")
