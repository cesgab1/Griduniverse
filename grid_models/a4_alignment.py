"""
A4 on the grid: does generic dynamics produce the alignment that gives tri-bimaximal neutrino mixing?
Two real A4 triplet 'flavon' fields phi (neutrino sector) and chi (charged-lepton sector).
Needed: phi along a body diagonal (1,1,1)-type, chi along an axis (1,0,0)-type (or vice versa; labels are arbitrary).
Method: build EVERY A4-invariant polynomial up to degree 4 in (phi, chi) with the Reynolds operator (group average),
give each an O(1) random coefficient (quadratic terms negative on average so the symmetry breaks), keep bounded potentials,
find the global minimum from many starts, and measure the misalignment angle of each field from its target direction.
Variant 'sequestered': terms that couple phi to chi suppressed by a factor s (0.1, 0.01): how much separation is needed.
"""
import numpy as np, sympy as sp, itertools, json, sys
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
S = np.diag([1, -1, -1]); T = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
G = [np.eye(3, dtype=int)]
while True:
    new = [g @ h for g in G for h in (S, T)]; added = False
    for n in new:
        if not any((n == g).all() for g in G): G.append(n); added = True
    if not added: break
print(f"A4 has {len(G)} elements")
x = sp.symbols("p1 p2 p3 c1 c2 c3"); phi, chi = sp.Matrix(x[:3]), sp.Matrix(x[3:])
inv = {2: [], 3: [], 4: []}
for deg in (2, 3, 4):
    seen = []
    for mono in itertools.combinations_with_replacement(x, deg):
        m = sp.Mul(*mono); acc = 0
        for g in G:
            gp, gc = sp.Matrix(g) * phi, sp.Matrix(g) * chi
            acc += m.subs(dict(zip(x, list(gp) + list(gc))), simultaneous=True)
        acc = sp.expand(acc / len(G))
        if acc == 0: continue
        P = sp.Poly(acc, *x); vec = dict(P.terms())
        seen.append((acc, vec))
    # linear independence
    keys = sorted({k for _, v in seen for k in v})
    Mx = np.array([[float(v.get(k, 0)) for k in keys] for _, v in seen])
    basis = []; B = np.zeros((0, len(keys)))
    for i, row in enumerate(Mx):
        if np.linalg.matrix_rank(np.vstack([B, row])) > B.shape[0]: B = np.vstack([B, row]); basis.append(seen[i][0])
    inv[deg] = basis
    print(f"degree {deg}: {len(basis)} independent A4 invariants")
def mixes(e):   # does the invariant couple phi and chi?
    fs = e.free_symbols; return any(s in fs for s in x[:3]) and any(s in fs for s in x[3:])
allinv = inv[2] + inv[3] + inv[4]
isq = np.array([d == 2 for d in [2]*len(inv[2]) + [3]*len(inv[3]) + [4]*len(inv[4])])
is4 = np.array([d == 4 for d in [2]*len(inv[2]) + [3]*len(inv[3]) + [4]*len(inv[4])])
mix = np.array([mixes(e) for e in allinv])
f_all = sp.lambdify([x], allinv, "numpy")
grad_all = sp.lambdify([x], [[sp.diff(e, s) for s in x] for e in allinv], "numpy")
def V(c, z): return float(np.dot(c, f_all(z)))
def dV(c, z): return np.array(grad_all(z), dtype=float).T @ c
DIAG = [np.array(v) / np.sqrt(3) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
AXES = [np.eye(3)[i] for i in range(3)]
def angle_to(v, targets):
    n = np.linalg.norm(v)
    if n < 1e-9: return np.pi / 2
    return min(np.arccos(min(1, abs(np.dot(v / n, t)))) for t in targets)
dirs = rng.normal(size=(3000, 6)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
def trial(s):
    c = rng.normal(0, 1, len(allinv))
    c[isq & ~mix] = -np.abs(rng.normal(1, 0.3, (isq & ~mix).sum()))        # negative mass terms for phi^2, chi^2
    c[mix] *= s
    q = np.array([np.dot(c * is4, f_all(d)) for d in dirs[:600]])           # quartic part must be positive in every direction
    if q.min() <= 0.02: return None
    best = None
    for _ in range(25):
        z0 = rng.normal(0, 1, 6)
        r = minimize(lambda z: V(c, z), z0, jac=lambda z: dV(c, z), method="BFGS", options=dict(gtol=1e-9))
        if best is None or r.fun < best.fun: best = r
    z = best.x; p, ch = z[:3], z[3:]
    if min(np.linalg.norm(p), np.linalg.norm(ch)) < 0.05 * max(np.linalg.norm(p), np.linalg.norm(ch)): return ("one field unbroken", None)
    m1 = max(angle_to(p, DIAG), angle_to(ch, AXES)); m2 = max(angle_to(ch, DIAG), angle_to(p, AXES))
    return ("ok", min(m1, m2))
out = {}
for s in [float(a) for a in (sys.argv[3].split(",") if len(sys.argv) > 3 else ["1.0", "0.1", "0.01"])]:
    mis = []; unb = 0; n = 0
    while n < int(sys.argv[2]) if len(sys.argv) > 2 else n < 300:
        t = trial(s)
        if t is None: continue
        n += 1
        if t[0] != "ok": unb += 1
        else: mis.append(t[1])
    mis = np.array(mis)
    fr = lambda a: np.mean(mis < a) * len(mis) / n
    print(f"cross-couplings x{s:<5}: {n} bounded random potentials | both fields broken in {len(mis)/n*100:.0f}% | "
          f"aligned within 0.02 rad: {fr(0.02)*100:.1f}%, within 0.15 rad: {fr(0.15)*100:.1f}%, within 0.3 rad: {fr(0.3)*100:.1f}% | median misalignment {np.median(mis):.2f} rad")
    out[s] = dict(n=n, broken_both=len(mis) / n, p002=fr(0.02), p015=fr(0.15), p03=fr(0.3), median=float(np.median(mis)), mis=mis.tolist())
TAG = sys.argv[3] if len(sys.argv) > 3 else 'all'
json.dump(out, open(f'/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/a4_alignment_{TAG}.json', 'w'))
