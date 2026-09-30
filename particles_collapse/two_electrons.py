"""
Two electrons on a grid (exact quantum calculations, 2D periodic grid L x L, hopping t = 1 between neighbours).
1. Same spin: antisymmetric two-particle state -> probability of finding them at separation d (the 'exchange hole').
2. Opposite spin with on-site cost U (Hubbard): chance both sit on the same point, vs U.
3. Grid that dimples (Holstein-type, adiabatic): each electron pulls its site down by u = g n / k; the dimple lowers the energy of
   any electron on that site by g u -> a shared dimple gives an on-site ATTRACTION  -lambda = -g^2/k  (net U_eff = U - lambda).
   Exact two-particle ground state vs lambda at fixed U: binding energy and pair size.
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, json
L = 12; N = L*L
def nbrs(i):
    x, y = divmod(i, L); return [((x+1)%L)*L+y, ((x-1)%L)*L+y, x*L+(y+1)%L, x*L+(y-1)%L]
T1 = sp.lil_matrix((N, N))
for i in range(N):
    for j in nbrs(i): T1[i, j] = -1.0
T1 = T1.tocsr(); I = sp.identity(N, format="csr")
H0 = sp.kron(T1, I) + sp.kron(I, T1)                       # two distinguishable particles on the grid
diag_same = np.zeros(N*N); diag_same[[i*N+i for i in range(N)]] = 1.0; D = sp.diags(diag_same)
def sep(i, j):
    xi, yi = divmod(i, L); xj, yj = divmod(j, L); dx = min(abs(xi-xj), L-abs(xi-xj)); dy = min(abs(yi-yj), L-abs(yi-yj)); return np.hypot(dx, dy)
S = np.array([[sep(i, j) for j in range(N)] for i in range(N)]).ravel()
out = {}
# --- 1. same spin: lowest antisymmetric state (project with swap operator)
P = sp.lil_matrix((N*N, N*N))
for i in range(N):
    for j in range(N): P[i*N+j, j*N+i] = 1.0
P = P.tocsr(); A = (sp.identity(N*N) - P)/2
# random-start power iteration in antisymmetric subspace via shifted Hamiltonian restricted
Ha = A @ H0 @ A + 20*(sp.identity(N*N) - A)
w, v = sla.eigsh(Ha, k=1, which="SA"); psi = v[:, 0]; prob = psi**2
dbins = [0, 1, np.sqrt(2), 2, 3, 4, 6]
print("1. SAME SPIN (antisymmetric): probability of being found at separation d, relative to independent particles")
rel_same = []
for d in dbins:
    m = np.isclose(S, d) if d < 4 else (S >= d)
    rel = prob[m].sum()/m.sum()*N*N; rel_same.append(rel)
    print(f"   d = {d:4.2f} grid spacings: {rel:.3f}")
out["same_spin"] = dict(d=dbins, rel=rel_same)
# --- 2. opposite spins, Hubbard U (symmetric spatial state)
print("\n2. OPPOSITE SPIN with on-site repulsion U: chance of sharing a grid point (independent particles: 1/N =", f"{1/N:.4f})")
res2 = []
for U in (0, 1, 2, 4, 8, 16, 1e3):
    w, v = sla.eigsh(H0 + U*D, k=1, which="SA"); p = v[:, 0]**2
    res2.append((U, float(p @ diag_same)))
    print(f"   U = {U:6g}: shared-point probability {p @ diag_same:.4f}  ({p @ diag_same*N:.2f} x independent)")
out["hubbard"] = res2
# --- 3. dimpling grid: net on-site interaction U - lambda
print("\n3. GRID THAT DIMPLES (electrons attract through the shared dimple), direct repulsion U = 4")
E1 = -4.0                                                   # single electron band bottom (2D, t=1)
res3 = []
for lam in (0, 2, 4, 5, 6, 8, 10, 14):
    Ueff = 4 - lam
    w, v = sla.eigsh(H0 + Ueff*D, k=1, which="SA"); p = v[:, 0]**2
    Eb = w[0] - 2*E1; size = np.sqrt(p @ S**2)
    res3.append((lam, float(Eb), float(size), float(p @ diag_same)))
    state = "BOUND PAIR" if Eb < -0.05 else "not bound (finite-grid level shift)" if Eb < 0 else "not bound"
    print(f"   dimple strength lambda = {lam:4.1f} (net {Ueff:+.1f}): pair energy vs two free electrons {Eb:+.3f} t, pair size {size:4.2f} spacings, on same point {p @ diag_same:.3f}  -> {state}")
out["dimple"] = res3
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/two_electrons.json", "w"), indent=1)
print(f"\n(free pair size on this {L}x{L} grid for comparison: {np.sqrt(np.mean(S**2)):.2f} spacings)")
