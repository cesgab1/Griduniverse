"""
ITERATION 42 (QG): does our random grid 'look different' (lower effective dimension) below the switch-over scale?
Effective (spectral) dimension from diffusion: a random walker's return probability P(sigma) after diffusion time sigma;
d_s = -2 dlnP/dln sigma. Spacetime walk on our grid with the preferred time slicing: time direction continuum (factor
sigma^-1/2), space = our random Mosaic (periodic box, points joined to neighbours, weights ~ 1/d^2, Laplacian eigenvalues lambda):
     P(sigma) ~ sigma^(-1/2) * mean_k exp(-sigma f(lambda_k)),   f = lambda (plain, z = 1)   or   f = lambda + lambda^3/M^4 (Horava z = 3)
Expectations written BEFORE running:
 D1 plain grid (z = 1): d_s ~ 4 at large scales, no plateau near 2; drops only at the single-cell scale (discreteness).
 D2 with the z = 3 rule switching on at a few cell spacings: a plateau d_s ~ 2.0 +/- 0.25 between the switch-over and the cell
    scale -- our random grid carries Horava's dimensional reduction cleanly (as CDT finds d_s ~ 2 at small scales).
 D3 (input, not tested) with z = 3 - eta, eta ~ 0.013: d_s = 1 + 3/z ~ 2.004 -- indistinguishable here.
The reduction itself is an input of the z = 3 rule (borrowed); what is tested is that our grid supports it over a real range.
"""
import numpy as np
from scipy.spatial import cKDTree
rng = np.random.default_rng(42)
n = 5000; L = n**(1/3)                          # mean spacing 1
X = rng.uniform(0, L, (n, 3)); tree = cKDTree(X, boxsize=L)
pairs = tree.query_pairs(1.6, output_type="ndarray")
d = X[pairs[:, 0]] - X[pairs[:, 1]]; d -= L*np.round(d/L); r = np.linalg.norm(d, axis=1)
W = np.zeros((n, n)); w = 1/r**2; W[pairs[:, 0], pairs[:, 1]] = w; W[pairs[:, 1], pairs[:, 0]] = w
Lap = np.diag(W.sum(1)) - W
lam = np.linalg.eigvalsh(Lap)
# calibrate to the continuum: the lowest non-zero modes should be (2 pi/L)^2 * |m|^2 with |m|^2 = 1 (6 modes)
cal = (2*np.pi/L)**2/np.median(lam[1:7]); lam = lam*cal
weyl = np.polyfit(np.log(lam[10:300]), np.log(np.arange(10, 300)), 1)[0]   # counting N(lambda) ~ lambda^(3/2) in 3-D
sig = np.geomspace(1e-3, 1e3, 400)
def ds(f):
    P = sig**-0.5*np.array([np.mean(np.exp(-s*f)) for s in sig]); return -2*np.gradient(np.log(P), np.log(sig))
out = ["ITERATION 42: spectral dimension of our random grid (expectations committed before running)", "",
       f"grid: {n} random points, periodic box, links within 1.6 spacings; Weyl-law exponent of the spectrum = {weyl:.2f} (3-D: 1.50)", ""]
for name, f in (("plain z = 1", lam), ("z = 3 switching on at ~3 spacings (M = 2)", lam + lam**3/2.0**4),
                ("z = 3 switching on at ~6 spacings (M = 1)", lam + lam**3/1.0**4)):
    D = ds(f); out.append(f"{name}:")
    for s_ in (300, 30, 3, 0.3, 0.03, 0.003):
        out.append(f"    diffusion scale sigma = {s_:7.3f} (~{np.sqrt(s_):5.2f} spacings): d_s = {np.interp(np.log(s_), np.log(sig), D):.2f}")
    mid = (sig > 0.01) & (sig < 1); out.append(f"    min / max d_s for sigma 0.01-1: {D[mid].min():.2f} / {D[mid].max():.2f}")
txt = "\n".join(out); print(txt); open("iter42_spectral_dimension.txt", "w").write(txt + "\n")
