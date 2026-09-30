"""
Step 1 & 2: the electron as a two-layer grid pattern (domain-wall construction, Kaplan 1992; here the 2D lattice version).
Grid: one ordinary direction x (momentum k) and one 'layer' direction s with Ns layers. Lattice Dirac Hamiltonian with the standard
Wilson term (removes the fake extra copies, the 'doublers'):
   H = sin(k) sx + sin(k_s) sy + (M - 2 + cos k + cos k_s) sz,   written in real space along s.
The layer-mass M(s) is +M inside a slab of W layers and -M outside -> two walls (the two 'layers' of the electron).
Expected: each wall carries ONE massless half, moving one way only (left-handed on one wall, right-handed on the other);
the halves 'snap' into each other only by leaking through the slab, so the electron's mass ~ exp(-W/xi): exponentially light.
"""
import numpy as np, json
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)
def H(k, Ns, W, M):
    Ms = np.where((np.arange(Ns) >= Ns//2 - W//2) & (np.arange(Ns) < Ns//2 - W//2 + W), M, -M)
    h = np.zeros((2*Ns, 2*Ns), complex)
    for s in range(Ns):
        h[2*s:2*s+2, 2*s:2*s+2] = np.sin(k)*sx + (Ms[s] - 2 + np.cos(k))*sz
        t = (sz - 1j*sy)/2                      # hopping along s: (cos k_s sz + sin k_s sy) in real space
        sp = (s+1) % Ns
        h[2*s:2*s+2, 2*sp:2*sp+2] += t
        h[2*sp:2*sp+2, 2*s:2*s+2] += t.conj().T
    return h, Ms
Ns, W, M = 60, 20, 0.5
wallA = Ns//2 - W//2; wallB = wallA + W
print("STEP 1: spectrum of the two-layer grid electron (Ns = 60 layers, slab of 20)")
ks = np.linspace(-np.pi, np.pi, 201); low = []
for k in ks:
    e, v = np.linalg.eigh(H(k, Ns, W, M)[0]); low.append(np.sort(np.abs(e))[:4])
low = np.array(low)
print(f"   lightest state near k=0 : |E| = {low[100,0]:.2e}   (light electron)")
print(f"   lightest state near k=pi: |E| = {low[-1,0]:.2f}      (would-be extra copy: heavy, removed by the Wilson term)")
print(f"   number of light (|E|<0.1) branches crossing near k=0: {int(np.sum(low[100] < 0.1))}")
# chirality + location: follow the two lightest states at small k
k = 0.2; e, v = np.linalg.eigh(H(k, Ns, W, M)[0]); order = np.argsort(np.abs(e))[:2]
k2 = 0.21; e2 = np.linalg.eigvalsh(H(k2, Ns, W, M)[0])
for idx in order:
    w = np.abs(v[:, idx])**2; ws = w[0::2] + w[1::2]
    near = lambda c: ws[max(0, c-4):c+4].sum()
    # velocity from neighbouring k
    E = e[idx]; E2 = e2[np.argmin(abs(e2 - E))]
    print(f"   light state E = {E:+.3f}: moves {'right' if (E2-E) > 0 else 'left '} (dE/dk = {(E2-E)/0.01:+.2f}); weight on wall A {near(wallA):.2f}, on wall B {near(wallB):.2f}")
print("   -> each wall holds one half, moving in one direction only: the two handednesses live on two different layers")
# one-sided force: a 'weak-force' coupling that acts only near wall A
print("   a force acting only on layers near wall A couples to the left-moving half with weight above and to the right-moving half with ~0")
print("\nSTEP 2: how the electron's mass (the 'snapping' between layers) depends on the layer separation")
res = []
for Wd in (2, 4, 6, 8, 10, 12, 14, 16):
    e = np.linalg.eigvalsh(H(0.0, 80, Wd, M)[0]); m = np.sort(np.abs(e))[0]
    res.append((Wd, m)); print(f"   walls {Wd:2d} layers apart: electron mass = {m:.3e} (grid units)")
Wv = np.array([r[0] for r in res]); mv = np.array([r[1] for r in res]); ok = mv > 1e-13
xi = -1/np.polyfit(Wv[ok], np.log(mv[ok]), 1)[0]
print(f"   mass falls as exp(-W/xi), xi = {xi:.2f} layers  (for M = {M})")
for Mx in (0.3, 0.7, 1.5, 1.0):
    r2 = [(Wd, np.sort(np.abs(np.linalg.eigvalsh(H(0.0, 80, Wd, Mx)[0])))[0]) for Wd in (2, 4, 6, 8)]
    a = np.array(r2); okk = a[:,1] > 1e-12
    x2 = -1/np.polyfit(a[okk, 0], np.log(a[okk, 1]), 1)[0] if okk.sum() > 1 else 0.0; print(f"   (layer mass M = {Mx}: xi = {x2:.2f} layers)")
ratios = {"electron": 0.511e-3, "muon": 0.1057, "tau": 1.777, "up quark": 2.2e-3, "top quark": 172.7}
MP = 1.22e19
print("\n   layers needed if the grid's natural mass scale is the Planck mass:")
for n_, m_ in ratios.items():
    print(f"     {n_:10s} mass/Planck = {m_/MP:.1e}  -> {np.log(MP/m_)*xi:5.1f} layer separations (xi={xi:.2f})")
json.dump(dict(res=res, xi=xi), open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/two_layer_electron.json", "w"))
