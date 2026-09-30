"""
Discrete version: each grid building block (a tetrahedron) carries a LABEL, not a continuous field:
   a in {0..3}: which vertex is occupied  -> contributes r_a   (diagonal direction)
   b in {0..5}: which edge is occupied    -> contributes m_b   (axis direction; opposite edges give opposite signs)
Blocks sit on a periodic 3D lattice (8^3). Random O(1) couplings, all local:
   J1 * [a_i == a_j]  (neighbours, same vertex)     J2 * [b_i == b_j]  (neighbours, same edge)
   K  * [a_i on b_i]  (same block, vertex on edge)   K2 * [a_i on b_j]  (neighbours)
Anneal with Metropolis to low temperature; the macroscopic flavons are phi = <r_a>, chi = <m_b>.
Measured: fraction of random coupling sets that end with phi exactly along a diagonal and chi exactly along an axis.
"""
import numpy as np, itertools, sys, json
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
R = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], float)
E = list(itertools.combinations(range(4), 2)); M = np.array([(R[a] + R[b]) / 2 for a, b in E])
ON = np.zeros((4, 6));
for i, (a, b) in enumerate(E): ON[a, i] = ON[b, i] = 1
L = 8; shifts = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
par = (np.indices((L, L, L)).sum(0) % 2)
DIAG = [r / np.sqrt(3) for r in R]; AXES = [np.eye(3)[i] for i in range(3)]
def ang(v, T):
    n = np.linalg.norm(v); return np.pi/2 if n < 1e-6 else min(np.arccos(min(1, abs(v @ t) / n)) for t in T)
def nb(A): return [np.roll(A, s, axis=(0, 1, 2)) for s in shifts]
def trial():
    J1, J2, K, K2 = rng.normal(0, 1, 4)
    a = rng.integers(0, 4, (L, L, L)); b = rng.integers(0, 6, (L, L, L))
    for T in np.geomspace(3.0, 0.05, 60):
        for _ in range(4):
            for p in (0, 1):
                m = par == p
                # vertex update
                na, nb_ = nb(a), nb(b)
                new = rng.integers(0, 4, (L, L, L))
                def Ea(x): return J1 * sum((x == q) for q in na) + K * ON[x, b] + K2 * sum(ON[x, q] for q in nb_)
                dE = Ea(new) - Ea(a); acc = m & (rng.random((L, L, L)) < np.exp(-np.clip(dE / T, -50, 50)))
                a = np.where(acc, new, a)
                na = nb(a)
                newb = rng.integers(0, 6, (L, L, L))
                def Eb(y): return J2 * sum((y == q) for q in nb_) + K * ON[a, y] + K2 * sum(ON[q, y] for q in na)
                dE = Eb(newb) - Eb(b); acc = m & (rng.random((L, L, L)) < np.exp(-np.clip(dE / T, -50, 50)))
                b = np.where(acc, newb, b)
    phi = R[a].mean(axis=(0, 1, 2)); chi = M[b].mean(axis=(0, 1, 2))
    return max(ang(phi, DIAG), ang(chi, AXES)), (J1, J2, K, K2), np.linalg.norm(phi), np.linalg.norm(chi)
n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
res = [trial() for _ in range(n)]
mis = np.array([r[0] for r in res]); J = np.array([r[1] for r in res])
print(f"discrete tetrahedral labels, {n} random coupling sets: aligned within 0.02 rad {np.mean(mis<0.02)*100:.0f}%, within 0.15 rad {np.mean(mis<0.15)*100:.0f}%, median {np.median(mis):.2f} rad")
ferro = (J[:, 0] < 0) & (J[:, 1] < 0)
print(f"   when both neighbour couplings favour equal labels ({ferro.mean()*100:.0f}% of sets): aligned within 0.15 rad {np.mean(mis[ferro] < 0.15)*100:.0f}%")
print(f"   otherwise: {np.mean(mis[~ferro] < 0.15)*100:.0f}%")
json.dump(dict(mis=mis.tolist(), J=J.tolist()), open(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/tetra_labels_{sys.argv[1] if len(sys.argv)>1 else 0}.json", "w"))
