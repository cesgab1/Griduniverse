"""
ITERATION 41b (declared follow-up to 41, same prediction): Ollivier-Ricci with COMMON RANDOM POINTS.
One flat random scatter (uniform in a ball); for the curved grid each point is moved radially so the scatter becomes uniform in
PROPER volume (proper volume inside the new radius = flat volume inside the old radius). The same point pairs are then measured
in both, so most sampling noise and bias cancel in Delta kappa = kappa_curved - kappa_flat.  Prediction (as in 41): eps^2 R(0)/30.
Expectations written BEFORE running:
 F1 if iteration 41's ~3x excess was sampling noise, the ratio here is within ~25% of 1 already at ~1000-2000 points per neighbourhood.
 F2 if it stays ~3, the excess is SYSTEMATIC (formula correction at delta = eps/2, or curvature varying across eps) and must be
    diagnosed, not attributed to noise.
"""
import numpy as np, ot
from scipy.spatial import cKDTree
from scipy.integrate import cumulative_trapezoid
A0, s = 0.3, 0.25; c0 = np.zeros(3)
def phi(x, A): return A*np.exp(-(x**2).sum(-1)/(2*s*s))
def plen(P, Q, A):
    t = np.linspace(0, 1, 5); w = np.array([1, 4, 2, 4, 1])/12.0
    return sum(wi*np.exp(phi(P + ti*(Q - P), A)) for ti, wi in zip(t, w))*np.linalg.norm(Q - P, axis=-1)
rg = np.linspace(0, 0.3, 30001)
Vc = cumulative_trapezoid(4*np.pi*rg**2*np.exp(3*A0*np.exp(-rg**2/(2*s*s))), rg, initial=0)
def to_curved(X):
    r = np.linalg.norm(X, axis=1); Vf = 4/3*np.pi*r**3; rn = np.interp(Vf, Vc, rg)
    return X*(rn/np.maximum(r, 1e-12))[:, None]
def kappa(X, A, i, j, tree, eps):
    def nb(k):
        c = np.array(tree.query_ball_point(X[k], eps*1.05)); d = plen(X[k][None], X[c], A); return c[(d <= eps) & (c != k)]
    a, b = nb(i), nb(j); delta = plen(X[i], X[j], A)
    C = plen(X[a][:, None, :], X[b][None, :, :], A)
    W = ot.emd2(np.full(len(a), 1/len(a)), np.full(len(b), 1/len(b)), C, numItermax=5000000)
    return 1 - W/delta, len(a)
eps = 0.07; R0 = 12*A0*np.exp(-2*A0)/s**2; pred = eps**2*R0/30
out = ["ITERATION 41b: Ollivier-Ricci with common random points (expectations committed before running)", "",
       f"prediction Delta kappa = {pred:.5f}", "", " points/neighbourhood   pairs   Delta kappa (mean +/- se)      ratio to prediction"]
rng = np.random.default_rng(411)
for rho, npairs in ((6e5, 60), (1.5e6, 40), (2.6e6, 24)):
    Rb = 0.2; n = rng.poisson(rho*4/3*np.pi*Rb**3); u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    Xf = u*(Rb*rng.uniform(0, 1, n)**(1/3))[:, None]; Xc = to_curved(Xf)
    tf, tc = cKDTree(Xf), cKDTree(Xc)
    centre = np.where(np.linalg.norm(Xf, axis=1) < 0.02)[0]; d = []; m = 0
    for i in rng.permutation(centre):
        c = np.array(tf.query_ball_point(Xf[i], 0.55*eps)); c = c[np.linalg.norm(Xf[c] - Xf[i], axis=1) > 0.45*eps]
        if len(c) == 0: continue
        j = c[0]
        kf, m = kappa(Xf, 0.0, i, j, tf, eps); kc, _ = kappa(Xc, A0, i, j, tc, eps); d.append(kc - kf)
        if len(d) >= npairs: break
    d = np.array(d); out.append(f"   ~{m:5d}            {len(d):3d}    {d.mean():+.5f} +/- {d.std()/np.sqrt(len(d)):.5f}          {d.mean()/pred:.2f} +/- {d.std()/np.sqrt(len(d))/pred:.2f}")
    print(out[-1], flush=True)
txt = "\n".join(out); open("iter41b_ollivier_same_points.txt", "w").write(txt + "\n")
