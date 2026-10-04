"""
ITERATION 43b (declared variance-reduction follow-up to 43; same criteria P1/P2): per grid point use SIX directions (+/- three
perpendicular axes, random orientation). The three axes sum to the full curvature R (Ric(e1)+Ric(e2)+Ric(e3) = R), and the
+/- pairs cancel the odd (gradient) terms that random single directions leave as scatter/bias where R changes quickly (edge).
"""
import numpy as np, ot, sys
from scipy.spatial.transform import Rotation
src = open("iter43_orc_profile.py").read().split("rf = np.linalg.norm(Xf, axis=1)")[0]
src = src.replace("rho = float(sys.argv[1]) if len(sys.argv) > 1 else 1.2e6", "rho = 1.2e6").replace(
      "npairs = int(sys.argv[2]) if len(sys.argv) > 2 else 30", "npairs = 0")
exec(src)
npts = int(sys.argv[1]) if len(sys.argv) > 1 else 8
rc = np.linalg.norm(Xc, axis=1)
for r0 in (0.0, 0.2, 0.35, 0.5):
    sel = np.where(np.abs(rc - r0) < 0.012)[0] if r0 > 0 else np.where(rc < 0.02)[0]
    est = []
    for i in rng.permutation(sel):
        Rm = Rotation.random(random_state=rng).as_matrix(); vals = []
        for e in np.r_[Rm, -Rm]:
            _, j = tf.query(Xf[i] + 0.5*eps*e)
            if j == i or not (0.4*eps < np.linalg.norm(Xf[j] - Xf[i]) < 0.6*eps): break
            kf, m = kap(Xf, 0.0, i, j, tf); kc, _ = kap(Xc, A0, i, j, tc); vals.append(kc - kf)
        if len(vals) == 6: est.append(30/eps**2*np.mean(vals))
        if len(est) >= npts: break
    est = np.array(est)
    line = f"r = {r0:4.2f}: measured R = {est.mean():+7.2f} +/- {est.std()/np.sqrt(len(est)):5.2f}   true R = {Rtrue(r0):+7.2f}   ({len(est)} points x 6 directions)"
    print(line, flush=True); open("iter43b_orc_triad.txt", "a").write(line + "\n")
