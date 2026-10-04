"""
ITERATION 41 (QG, Coalesce's question): is the tidy-cell requirement only an artefact of Regge's curvature DEFINITION?
Test a definition built for random grids: Ollivier-Ricci curvature (ORC). For two grid points x, y a distance delta apart,
compare the 'cost' W of moving the cloud of x's neighbours (within eps) onto y's neighbours:
    kappa(x, y) = 1 - W / delta.     Curved space shifts kappa by  eps^2 Ric(v,v) / (2 (n + 2))  (Ollivier 2009; n = 3: /10).
Random grid: points sprinkled uniformly in PROPER volume (density ~ sqrt(g)), lengths of all short links from the metric
g = e^{2 phi} delta with a bump phi = A exp(-r^2/2s^2) (A = 0.3, s = 0.25). At the centre, by symmetry, Ric(v,v) = R/3 with
R(0) = 12 A e^{-2A} / s^2.  Measurement: mean kappa over many pairs near the centre, curved MINUS flat (same procedure, A = 0),
so the finite-sampling bias of ORC cancels.  Prediction: Delta kappa = eps^2 R(0) / 30.
Expectations written BEFORE running:
 O1 Delta kappa approaches the prediction as the number of grid points per neighbourhood grows; within ~25% at the densest run.
 O2 its scatter shrinks with density.
 Kill: no approach to the prediction (wrong sign, or off by > 50% at the densest run) -> the averaging definition does not
 rescue the random grid either.
"""
import numpy as np, ot
from scipy.spatial import cKDTree
A0, s = 0.3, 0.25; c0 = np.array([0.5, 0.5, 0.5])
def phi(x, A): return A*np.exp(-((x - c0)**2).sum(-1)/(2*s*s))
def plen(P, Q, A):
    t = np.linspace(0, 1, 5); w = np.array([1, 4, 2, 4, 1])/12.0
    return sum(wi*np.exp(phi(P + ti*(Q - P), A)) for ti, wi in zip(t, w))*np.linalg.norm(Q - P, axis=-1)
def sprinkle(rho, A, rng, half=0.17):
    V = (2*half)**3; n = rng.poisson(rho*V*np.exp(3*A))
    X = c0 + rng.uniform(-half, half, (n, 3)); keep = rng.uniform(0, 1, n) < np.exp(3*phi(X, A))/np.exp(3*A)
    return X[keep]
def mean_kappa(rho, A, eps, rng, npairs=120):
    X = sprinkle(rho, A, rng); tree = cKDTree(X)
    centre = np.where(np.linalg.norm(X - c0, axis=1) < 0.02)[0]; ks = []
    for i in rng.permutation(centre):
        cand = tree.query_ball_point(X[i], eps*1.05)      # coordinate radius >= proper radius (phi >= 0)
        cand = np.array(cand); d = plen(X[i][None], X[cand], A); nb_x = cand[(d <= eps) & (cand != i)]
        dy = plen(X[i][None], X[nb_x], A); ys = nb_x[(dy > 0.4*eps) & (dy < 0.6*eps)]
        if len(ys) == 0: continue
        j = ys[0]; delta = plen(X[i], X[j], A)
        cj = np.array(tree.query_ball_point(X[j], eps*1.05)); nb_y = cj[(plen(X[j][None], X[cj], A) <= eps) & (cj != j)]
        if len(nb_x) < 5 or len(nb_y) < 5: continue
        C = plen(X[nb_x][:, None, :], X[nb_y][None, :, :], A)
        W = ot.emd2(np.full(len(nb_x), 1/len(nb_x)), np.full(len(nb_y), 1/len(nb_y)), C)
        ks.append(1 - W/delta)
        if len(ks) >= npairs: break
    return np.mean(ks), np.std(ks)/np.sqrt(len(ks)), len(ks), len(nb_x)
eps = 0.07
R0 = 12*A0*np.exp(-2*A0)/s**2; pred = eps**2*R0/30
out = ["ITERATION 41: Ollivier-Ricci curvature on the RANDOM grid (expectations committed before running)", "",
       f"R(0) = {R0:.2f};  prediction Delta kappa = eps^2 R/30 = {pred:.5f}  (eps = {eps})", "",
       " points/neighbourhood   kappa curved          kappa flat            Delta kappa (curved - flat)   ratio to prediction"]
rng = np.random.default_rng(41)
for rho in (2e5, 6e5, 1.5e6):
    kc, ec, nc, mc = mean_kappa(rho, A0, eps, rng); kf, ef, nf, mf = mean_kappa(rho, 0.0, eps, rng)
    dk = kc - kf; err = np.hypot(ec, ef)
    out.append(f"   ~{mc:5d}            {kc:+.5f}+/-{ec:.5f}   {kf:+.5f}+/-{ef:.5f}   {dk:+.5f} +/- {err:.5f}           {dk/pred:.2f} +/- {err/pred:.2f}")
txt = "\n".join(out); print(txt); open("iter41_ollivier_ricci.txt", "w").write(txt + "\n")
