"""
Is a light electron natural on the grid? Shake piece 16's link phases randomly (a stand-in for quantum fluctuations of the field, flux Q = 0)
and measure the mass the fermion picks up with no mass put in.  Mass = smallest |Re(eigenvalue)| of the low modes after mapping onto the
continuum axis (overlap: stereographic map of the Ginsparg-Wilson circle).  Chiral symmetry forbids any induced mass.
"""
import numpy as np
exec(open("fermions_index.py").read().split("L, r, m0")[0])
rng2 = np.random.default_rng(4)
L, r, m0 = 20, 1.0, 1.0
n, P, links, Vc = grid(L); G5 = np.kron(np.eye(n), sz)
Dn0, W0, _, _ = operators(n, P, links, Vc, L, 0)
pair = {}
print("phase noise sd | Wilson: induced mass of lowest mode | overlap: induced mass | GW residual")
for sd in [0.0, 0.3, 0.6, 0.9]:
    th = {}
    for a, b, A, d, s in links:
        key = (min(a, b), max(a, b), round(abs(d[0]), 6), round(abs(d[1]), 6))
        if key not in th: th[key] = rng2.normal(0, sd)
    ph = np.ones((n, n), complex); Kmask = np.zeros((n, n))
    Dn = Dn0.copy(); W = W0.copy()
    for a, b, A, d, s in links:
        key = (min(a, b), max(a, b), round(abs(d[0]), 6), round(abs(d[1]), 6)); t = th[key] * (1 if a < b else -1)
        f = np.exp(1j * t)
        Dn[2*a:2*a+2, 2*b:2*b+2] *= f; W[2*a:2*a+2, 2*b:2*b+2] *= f
    DW = Dn + r / 2 * W
    ew = np.linalg.eigvals(DW); low = np.argsort(np.abs(ew))[:4]
    H = G5 @ (DW - m0 * np.eye(2 * n)); H = (H + H.conj().T) / 2
    e, U = eigh(H); Dov = m0 * (np.eye(2 * n) + G5 @ ((U * np.sign(e)) @ U.conj().T))
    gw = np.linalg.norm(G5 @ Dov + Dov @ G5 - Dov @ G5 @ Dov / m0) / np.linalg.norm(Dov)
    eo = np.linalg.eigvals(Dov); eo = eo[np.abs(eo - 2 * m0) > 1e-6]; ep = eo / (1 - eo / (2 * m0)); lo = np.argsort(np.abs(ep))[:4]
    print(f"   {sd:4.1f}        |   {np.abs(ew[low].real).min():.4f}  (lowest |E| {np.abs(ew[low]).min():.4f})      |   {np.abs(ep[lo].real).max():.1e}  (lowest |E| {np.abs(ep[lo]).min():.4f}) | {gw:.0e}")
