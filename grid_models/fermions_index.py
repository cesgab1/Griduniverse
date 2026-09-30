"""
Spin-1/2 + electromagnetism on the random grid: the index theorem.
Put Q whole units of magnetic flux through the periodic random grid using piece 16's link phases,
then count zero-energy states of the overlap operator and their handedness (gamma5 = sigma_z).
Continuum (Atiyah-Singer): exactly |Q| zero modes, all of one handedness; no zero modes for Q = 0.
"""
import numpy as np, json
from scipy.spatial import Voronoi
from scipy.linalg import eigh
rng = np.random.default_rng(12)
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)

def poly_area(V):
    c = V.mean(0); a = np.arctan2(V[:, 1] - c[1], V[:, 0] - c[0]); V = V[np.argsort(a)]
    return 0.5 * abs(np.dot(V[:, 0], np.roll(V[:, 1], -1)) - np.dot(V[:, 1], np.roll(V[:, 0], -1)))

def grid(L):
    n = rng.poisson(L * L); P = rng.uniform(0, L, (n, 2))
    shifts = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)])
    Q = np.concatenate([P + s * L for s in shifts]); home = np.tile(np.arange(n), 9); sh = np.repeat(shifts, n, axis=0); c0 = 4 * n
    vor = Voronoi(Q); links = []
    for (a, b), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in rv: continue
        A = np.linalg.norm(vor.vertices[rv[0]] - vor.vertices[rv[1]])
        for s, t in ((a, b), (b, a)):
            if c0 <= s < c0 + n: links.append((home[s], home[t], A, Q[t] - Q[s], sh[t]))
    links = list({(a, b, round(d[0], 8), round(d[1], 8)): (a, b, A, d, s) for a, b, A, d, s in links}.values())
    Vc = np.array([poly_area(vor.vertices[vor.regions[vor.point_region[c0 + i]]]) for i in range(n)])
    return n, P, links, Vc

def operators(n, P, links, Vc, L, Qtop):
    B = 2 * np.pi * Qtop / L**2
    K = np.zeros((2 * n, 2 * n), complex); Lp = np.zeros((n, n), complex)
    for a, b, A, d, s in links:
        xa, ya = P[a]; xb = P[b][0]
        # transporter from the ghost copy of b back to a, A = (-B y, 0) on the plane, plus the torus patching e^{-i s_y B L x_b}
        U = np.exp(1j * B * (ya + d[1] / 2) * d[0]) * np.exp(-1j * s[1] * B * L * xb)
        nv = d / np.linalg.norm(d)
        K[2*a:2*a+2, 2*b:2*b+2] += (A / 2) * U * (nv[0] * sx + nv[1] * sy)
        w = A / np.linalg.norm(d); Lp[a, b] -= w * U; Lp[a, a] += w
    herm_err = np.linalg.norm(K + K.conj().T) / np.linalg.norm(K)       # must be ~0 if the phases are consistent
    s_ = np.repeat(1 / np.sqrt(Vc), 2)
    Dn = s_[:, None] * K * s_[None, :]
    W = np.kron(Lp / np.sqrt(np.outer(Vc, Vc)), np.eye(2))
    return Dn, W, herm_err, np.linalg.norm(Lp - Lp.conj().T) / np.linalg.norm(Lp)

L, r, m0 = 24, 1.0, 1.0
out = []
for trial in range(2):
    n, P, links, Vc = grid(L)
    G5 = np.kron(np.eye(n), sz)
    for Qtop in [0, 1, 2, 3, -2]:
        Dn, W, e1, e2 = operators(n, P, links, Vc, L, Qtop)
        H = G5 @ (Dn + r / 2 * W - m0 * np.eye(2 * n)); H = (H + H.conj().T) / 2
        e, U = eigh(H); Dov = m0 * (np.eye(2 * n) + G5 @ ((U * np.sign(e)) @ U.conj().T))
        HO = Dov.conj().T @ Dov; HO = (HO + HO.conj().T) / 2
        lam, V = eigh(HO)
        zero = lam < 1e-8
        chir = np.real(np.einsum("ij,ik,kj->j", V[:, zero].conj(), G5, V[:, zero]))
        index = -0.5 * np.trace(G5 @ Dov).real / m0                # lattice topological charge from the operator itself
        nxt = np.sqrt(lam[~zero][0])
        print(f"grid {trial} ({n} pts)  flux Q = {Qtop:+d}: zero modes {zero.sum()}, handedness {np.round(chir, 3).tolist()}, "
              f"index from trace {index:+.4f}, next level {nxt:.3f}  [phase consistency {e1:.0e}, {e2:.0e}]")
        out.append(dict(trial=trial, Q=Qtop, zeros=int(zero.sum()), chir=chir.tolist(), index=float(index)))
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/fermion_index_results.json", "w"), indent=1)
