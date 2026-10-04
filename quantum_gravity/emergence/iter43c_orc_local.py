"""
ITERATION 43c (declared fix of 43/43b; same criteria P1/P2). Diagnosis of the edge failure (r = 0.5: -50 and -35 vs true -2.6):
the GLOBAL volume-matching map moved points by ~0.08 and sheared the pattern by ~40% at r = 0.5 (all of the bump's extra proper
volume pushes the outer points inward), while the true curvature there is tiny -> the shear, not curvature, dominated Delta kappa.
Fix: a LOCAL map around each measurement point x0: Y = x0 + (X - x0) e^{-phi(x0)} (isotropic rescale, no shear). Proper density is
then uniform near x0 to first order; the leftover gradient error is odd in direction and cancels in the +/- six-direction average.
"""
import numpy as np, ot, sys
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
A0, s = 0.3, 0.25
def phi(x, A): return A*np.exp(-(x**2).sum(-1)/(2*s*s))
def plen(P, Q, A):
    t = np.linspace(0, 1, 5); w = np.array([1, 4, 2, 4, 1])/12.0
    return sum(wi*np.exp(phi(P + ti*(Q - P), A)) for ti, wi in zip(t, w))*np.linalg.norm(Q - P, axis=-1)
def Rtrue(r):
    ph = A0*np.exp(-r**2/(2*s*s)); lap = ph*(r**2/s**4 - 3/s**2); g2 = ph**2*r**2/s**4
    return -np.exp(-2*ph)*(4*lap + 2*g2)
eps = 0.07; rho = 1.2e6; rng = np.random.default_rng(431)
shells = [float(a) for a in sys.argv[1].split(",")] if len(sys.argv) > 1 else [0.35, 0.5]
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for r0 in shells:
    est = []
    while len(est) < npts:
        u = rng.normal(size=3); x0 = r0*u/np.linalg.norm(u) if r0 > 0 else np.zeros(3)
        box = 0.2; m_ = rng.poisson(rho*(2*box)**3); Xf = x0 + rng.uniform(-box, box, (m_, 3)); Xf[0] = x0
        sc = np.exp(-phi(x0, A0)); Y = x0 + (Xf - x0)*sc
        tf, tc = cKDTree(Xf), cKDTree(Y)
        def kap(X, tree, A, i, j, scale):
            def nb(k):
                c = np.array(tree.query_ball_point(X[k], eps*1.3*scale + 1e-9)); d = plen(X[k][None], X[c], A); return c[(d <= eps) & (c != k)]
            a, b = nb(i), nb(j); C = plen(X[a][:, None, :], X[b][None, :, :], A)
            return 1 - ot.emd2(np.full(len(a), 1/len(a)), np.full(len(b), 1/len(b)), C, numItermax=5000000)/plen(X[i], X[j], A)
        Rm = Rotation.random(random_state=rng).as_matrix(); vals = []
        for e in np.r_[Rm, -Rm]:
            _, j = tf.query(x0 + 0.5*eps*e)
            if j == 0 or not (0.4*eps < np.linalg.norm(Xf[j] - x0) < 0.6*eps): break
            vals.append(kap(Y, tc, A0, 0, j, sc) - kap(Xf, tf, 0.0, 0, j, 1.0))
        if len(vals) == 6: est.append(30/eps**2*np.mean(vals))
    est = np.array(est)
    line = f"r = {r0:4.2f}: measured R = {est.mean():+7.2f} +/- {est.std()/np.sqrt(len(est)):5.2f}   true R = {Rtrue(r0):+7.2f}   ({len(est)} points x 6 directions, local map)"
    print(line, flush=True); open("iter43c_orc_local.txt", "a").write(line + "\n")
