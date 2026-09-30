"""
Graded drainage: a fraction D of links is drained (loose, carries nothing); the rest are sealed (normal springs).
Question 1 (pure geometry, no inputs): how does the patch stiffness K depend on the sealed fraction f?
   Random Delaunay mosaics, 2D and 3D, periodic, linear response (bulk + shear), several random draws.
Question 2: if drainage falls linearly with pull (uniform spread of throat strengths), f = f0 + (1-f0) x,
   what law g(g_N) results, and does it give the square-root law?   sigma(x) = Int_0^x K(f(x')) dx'.
"""
import numpy as np, json
from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from itertools import product

def build(N, d, rng):
    u = rng.random((N, d)); shifts = np.array(list(product([-1, 0, 1], repeat=d)))
    P = np.concatenate([u + s for s in shifts]); idx = np.tile(np.arange(N), len(shifts)); sh = np.repeat(shifts, N, axis=0)
    E = {}
    for simp in Delaunay(P).simplices:
        for a in range(d+1):
            for b in range(a+1, d+1):
                i, j = simp[a], simp[b]
                if not (np.all(sh[i] == 0) or np.all(sh[j] == 0)): continue
                ii, jj, s = idx[i], idx[j], sh[j] - sh[i]
                if (ii, tuple(s)) > (jj, tuple(-s)): ii, jj, s = jj, ii, -s
                E[(ii, jj, tuple(s))] = 1
    k = np.array(list(E.keys()), dtype=object)
    I = np.array([e[0] for e in E]); J = np.array([e[1] for e in E]); S = np.array([e[2] for e in E], float)
    dv = u[J] - u[I] + S; L = np.linalg.norm(dv, axis=1)
    return u, I, J, dv/L[:, None], L

def moduli(N, d, I, J, n, L, keep):
    I, J, n, L = I[keep], J[keep], n[keep], L[keep]; kk = 1/L         # EA = 1
    rows, cols, vals = [], [], []
    for a in range(d):
        for b in range(d):
            v = kk*n[:, a]*n[:, b]
            for (p, q, sgn) in ((I, I, 1), (J, J, 1), (I, J, -1), (J, I, -1)):
                rows.append(p*d+a); cols.append(q*d+b); vals.append(sgn*v)
    H = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N*d, N*d)).tocsr()
    H = H + 1e-9*np.mean(kk)*__import__("scipy.sparse", fromlist=["eye"]).eye(N*d)
    out = []
    strains = [np.eye(d)/d]                                             # bulk
    Sh = np.zeros((d, d)); Sh[0, 1] = Sh[1, 0] = 0.5; strains.append(Sh) # shear
    for Eps in strains:
        ea = np.einsum("ia,ab,ib->i", n, Eps, n)*L                      # affine extension
        fa = kk*ea                                                      # affine tension
        F = np.zeros((N, d)); np.add.at(F, I, (fa[:, None]*n)); np.add.at(F, J, -(fa[:, None]*n))
        du = spsolve(H, F.ravel()).reshape(N, d)                       # nonaffine relaxation
        ext = ea + np.einsum("ia,ia->i", du[J] - du[I], n)
        out.append(np.sum(kk*ext**2))                                   # 2 * energy (volume = 1)
    return out

res = {}
for d, N, seeds in ((2, 1500, 4), (3, 600, 3)):
    fs = np.round(np.concatenate([np.arange(0.40, 0.95, 0.025), [0.95, 0.975, 1.0]]), 3)
    Kb = np.zeros((seeds, len(fs))); Ks = np.zeros_like(Kb)
    for sd in range(seeds):
        rng = np.random.default_rng(100+sd); u, I, J, n, L = build(N, d, rng)
        order = rng.permutation(len(L))
        for m, f in enumerate(fs):
            keep = np.zeros(len(L), bool); keep[order[:int(round(f*len(L)))]] = True
            Kb[sd, m], Ks[sd, m] = moduli(N, d, I, J, n, L, keep)
    Kb /= Kb[:, -1:]; Ks /= Ks[:, -1:]
    z = 2*len(L)/N
    res[d] = dict(f=fs.tolist(), bulk=Kb.mean(0).tolist(), shear=Ks.mean(0).tolist(), z_full=z)
    print(f"\n{d}D Delaunay mosaic, full coordination z = {z:.2f}; isostatic z = {2*d} -> naive rigidity loss at f = {2*d/z:.3f}")
    print("   sealed f   bulk K/K_full   shear G/G_full")
    for m, f in enumerate(fs): print(f"   {f:6.3f}     {Kb[:, m].mean():.4f}         {Ks[:, m].mean():.4f}")
json.dump(res, open("drained_stiffness.json", "w"))
