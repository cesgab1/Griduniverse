"""
ITERATION 44 (QG): the dark-energy law's ingredients from fields living ON the random grid (not the continuum toy).
Random Mosaic (periodic box, Laplacian as in iteration 42). Two lattice fields:
   Phi : the tension field (minimally coupled, massless); tension = its velocity Pi.  Pi' = -3H Pi + Lap(Phi) + g (chi^2 - <chi^2>)
   chi : the faint thermal relic (its own lattice wave equation, started in a thermal state), which jostles Phi through g chi^2.
Background: constant H (slow compared with lattice frequencies). Measured:
   M  memory: autocorrelation time of the UNIFORM part of Pi (box average) -> expected 1/(3H)  (Hubble friction)
   V  variance of Pi averaged over regions of n cells -> expected ~ 1/n for regions larger than the relic's correlation length
Expectations written BEFORE running: M within 15% of 1/(3H); V slope -1.0 +/- 0.15 over the larger regions.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix, diags
rng = np.random.default_rng(44)
n = 3000; L = n**(1/3); X = rng.uniform(0, L, (n, 3)); tree = cKDTree(X, boxsize=L)
pr = tree.query_pairs(1.6, output_type="ndarray"); d = X[pr[:, 0]] - X[pr[:, 1]]; d -= L*np.round(d/L); w = 1/(d**2).sum(1)
W = csr_matrix((np.r_[w, w], (np.r_[pr[:, 0], pr[:, 1]], np.r_[pr[:, 1], pr[:, 0]])), shape=(n, n))
Lap = diags(np.asarray(W.sum(1)).ravel()) - W                    # positive semi-definite
H = 0.01; g = 0.05; dt = 0.05; T = 1.0
# thermal relic: random positions/velocities with equipartition (classical thermal state of the lattice field)
chi = rng.normal(0, np.sqrt(T), n); chiv = rng.normal(0, np.sqrt(T), n)
Phi = np.zeros(n); Pi = np.zeros(n)
steps = 120000; rec = 20; Pis = []
for k in range(steps):
    src = g*(chi**2 - np.mean(chi**2))
    Pi += dt*(-3*H*Pi - Lap @ Phi + src); Phi += dt*Pi
    chiv += dt*(-(Lap @ chi) - 0.0*chi); chi += dt*chiv
    if k % rec == 0: Pis.append(Pi.copy())
P = np.array(Pis)[len(Pis)//5:]; tstep = dt*rec
u = P.mean(1); u -= u.mean(); ac = np.correlate(u, u, "full")[len(u) - 1:]; ac /= ac[0]
tm = np.argmax(ac < np.exp(-1))*tstep
# region averages: spheres of growing radius around random centres
cent = X[rng.choice(n, 40, replace=False)]; rows = []
for R in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
    vs, ns = [], []
    for c in cent:
        dd = X - c; dd -= L*np.round(dd/L); idx = np.where((dd**2).sum(1) < R*R)[0]
        if len(idx) < 3: continue
        reg = P[:, idx].mean(1); vs.append(np.var(reg - u - u.mean() + 0*reg)); ns.append(len(idx))
    rows.append((np.mean(ns), np.mean(vs)))
rows = np.array(rows); slope = np.polyfit(np.log(rows[-4:, 0]), np.log(rows[-4:, 1]), 1)[0]
out = ["ITERATION 44: tension field jostled by a thermal relic, both living on the random grid (expectations committed first)", "",
       f"memory of the uniform tension: {tm:.1f} time units; Hubble-friction expectation 1/(3H) = {1/(3*H):.1f}  (ratio {tm*3*H:.2f})",
       " region size (cells)   variance of region-averaged tension (uniform part removed)"]
for nn, vv in rows: out.append(f"   {nn:7.1f}            {vv:.3e}")
out.append(f"slope over the larger regions: {slope:+.2f} (expected -1.0: independent cells -> sqrt(N) averaging)")
txt = "\n".join(out); print(txt); open("iter44_lattice_jostle.txt", "w").write(txt + "\n")
