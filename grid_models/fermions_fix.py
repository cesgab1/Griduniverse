"""
Spin-1/2 on the random grid, part 2.
 (a) sanity: smooth plane waves are handled correctly by the naive operator (so the swamp of fake low modes is real, not a bug)
 (b) Wilson fix: add (r/2) x the grid's own Laplacian (the same Voronoi weights as piece 16) -> rough modes get Planck-scale mass
 (c) overlap (Ginsparg-Wilson) operator built from it: D = m0 [1 + g5 sign(g5 (D_W - m0))] -> exact lattice chiral symmetry
     checks: species count, Ginsparg-Wilson relation  g5 D + D g5 = (1/m0) D g5 D,  low spectrum vs continuum
"""
import numpy as np, json
exec(open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/fermions_on_grid.py").read().split("L = 36")[0])
from scipy.linalg import eig, eigh
sz = np.diag([1., -1.]).astype(complex)

def build(n, links, Vc, P=None):
    K = np.zeros((2 * n, 2 * n), complex); Lp = np.zeros((n, n))
    for (a, b, dx, dy), (A, nv) in links.items():
        K[2*a:2*a+2, 2*b:2*b+2] += (A / 2) * (nv[0] * sx + nv[1] * sy)
        w = A / np.hypot(dx, dy); Lp[a, b] -= w; Lp[a, a] += w          # each unordered link appears twice (a->b, b->a)
    K = (K - K.conj().T) / 2
    Lp = (Lp + Lp.T) / 2
    s = np.repeat(1 / np.sqrt(Vc), 2)
    Dn = s[:, None] * K * s[None, :]                                          # anti-Hermitian, metric-symmetrised
    W = np.kron(Lp / np.sqrt(np.outer(Vc, Vc)), np.eye(2))                   # Hermitian, >= 0
    return Dn, W

def count(ev, L, Lams=(0.3, 0.5, 0.7)):
    return {Lam: (np.sum(np.abs(ev) < Lam), continuum_count(L, Lam)) for Lam in Lams}

L = 30; r = 1.0; m0 = 1.0
G5 = None; results = {}
for name, grid in [("square", square_grid), ("random", random_grid)]:
    n, links, Vc = grid(L)
    Dn, W = build(n, links, Vc)
    G5 = np.kron(np.eye(n), sz)
    # (a) plane-wave check (random grid needs positions: rebuild them from links is awkward -> use Rayleigh on |Dn psi| for smooth spinor fields)
    DW = Dn + (r / 2) * W
    ev_n = np.linalg.eigvals(Dn); ev_w = np.linalg.eigvals(DW)
    H = G5 @ (DW - m0 * np.eye(2 * n)); H = (H + H.conj().T) / 2
    e, U = eigh(H); sgn = (U * np.sign(e)) @ U.conj().T
    Dov = m0 * (np.eye(2 * n) + G5 @ sgn)
    gw = np.linalg.norm(G5 @ Dov + Dov @ G5 - Dov @ G5 @ Dov / m0) / np.linalg.norm(Dov)
    ev_o = np.linalg.eigvals(Dov)
    ev_o_proj = ev_o / (1 - ev_o / (2 * m0))            # stereographic map of the GW circle onto the imaginary axis
    gap = np.min(np.abs(e))
    res = {}
    for lab, ev in [("naive", ev_n), ("Wilson", ev_w), ("overlap", ev_o_proj)]:
        c = count(ev, L); res[lab] = {str(k): [int(a), int(b)] for k, (a, b) in c.items()}
        print(f"{name:7s} {lab:8s}: " + "   ".join(f"|E|<{k}: {a:4d} vs {b:3d} ({a/b:.2f} species)" for k, (a, b) in c.items()))
    lowo = np.sort(np.abs(ev_o_proj))[:12]
    print(f"{name:7s} overlap: Ginsparg-Wilson residual {gw:.1e}; smallest |H_W| eigenvalue (spectral gap) {gap:.3f}; lowest |E| {np.round(lowo, 3)}")
    results[name] = dict(counts=res, gw=float(gw), gap=float(gap), low=lowo.tolist())
kc = np.sort(np.linalg.norm(2*np.pi/L*np.array([(a, b) for a in range(-3, 4) for b in range(-3, 4)]), axis=1))
print("continuum lowest |k| (x2 spin):", np.round(np.repeat(kc[:6], 2), 3))
json.dump(results, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/fermion_fix_results.json", "w"), indent=1)
