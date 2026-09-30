"""
Probability version of the lepton census (the standard way to judge 'anarchic' neutrino textures):
for every grid-compatible (eps, h, L) with at least one singlet choice e that fits e, mu, tau within x3 (median),
draw 400 sets of O(1) coefficients and count how often ALL of these land inside the measured ranges:
  tight (~3 sigma): sin^2 th12 in [0.27,0.35], sin^2 th23 in [0.41,0.62], sin^2 th13 in [0.019,0.025], dm21/dm31 in [0.024,0.036]
  loose (x2):       [0.2,0.45], [0.35,0.65], [0.01,0.05], [0.01,0.08]
Baseline for comparison: pure anarchy (all lepton doublet charges equal).
"""
import numpy as np, itertools, glob, json, sys
from math import gcd
rng = np.random.default_rng(7)
v = 246.22; y_obs = np.log(np.array([0.4866e-3, 0.1027, 1.746]) * np.sqrt(2) / v)
def zmin(n, N): return np.minimum(n % N, N - n % N)
ND0, ND = 24, 400
ce0 = np.exp(rng.uniform(np.log(.5), np.log(2), (ND0, 3, 3))) * np.exp(2j*np.pi*rng.random((ND0, 3, 3)))
ce = np.exp(rng.uniform(np.log(.5), np.log(2), (ND, 3, 3))) * np.exp(2j*np.pi*rng.random((ND, 3, 3)))
cn = np.exp(rng.uniform(np.log(.5), np.log(2), (ND, 3, 3))) * np.exp(2j*np.pi*rng.random((ND, 3, 3))); cn = (cn + np.transpose(cn, (0, 2, 1))) / 2
dm31, dm21 = 2.51e-3, 7.42e-5
def neutrino_stats(nn, ne, eps):
    Un, sn, _ = np.linalg.svd(cn * eps**nn); Un = Un[:, :, ::-1]; sn = sn[:, ::-1]
    Ue, _, _ = np.linalg.svd(ce * eps**ne); Ue = Ue[:, :, ::-1]
    P = np.abs(np.einsum("dji,djl->dil", Ue.conj(), Un))**2
    s13 = P[:, 0, 2]; s12 = P[:, 0, 1] / (1 - s13 + 1e-12); s23 = P[:, 1, 2] / (1 - s13 + 1e-12)
    r = (sn[:, 1]**2 - sn[:, 0]**2) / (sn[:, 2]**2 - sn[:, 0]**2 + 1e-300)
    tight = (s12 > .27) & (s12 < .35) & (s23 > .41) & (s23 < .62) & (s13 > .019) & (s13 < .025) & (r > .024) & (r < .036)
    loose = (s12 > .2) & (s12 < .45) & (s23 > .35) & (s23 < .65) & (s13 > .01) & (s13 < .05) & (r > .01) & (r < .08)
    ratio = sn[:, 0] / sn[:, 2]; m3 = np.sqrt(dm31 / (1 - ratio**2)); msum = ratio*m3 + np.sqrt((ratio*m3)**2 + dm21) + m3
    return tight.mean(), loose.mean(), msum[loose] if loose.any() else np.array([])
# baseline: anarchy
nn = np.zeros((3, 3), int); t, l, ms = neutrino_stats(nn, np.array([[4, 2, 0]] * 3).T * 0 + np.array([[5, 2, 0]]), 0.1)
print(f"baseline, pure anarchy (equal lepton-doublet charges): tight {t*100:.2f}%, loose {l*100:.1f}%, median sum m_nu of loose hits {np.median(ms) if len(ms) else float('nan'):.3f} eV")
Ns = [int(a) for a in sys.argv[1].split(",")]
for f in sorted(glob.glob("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/flav_ZN_*_nu.json")):
    N = int(f.split("_")[-2])
    if N not in Ns: continue
    Q = np.array(json.load(open(f))["charges"][0:3]); g = gcd(9, N); need = (-3 * Q.sum()) % g
    Ls = [np.array(L) for L in itertools.combinations_with_replacement(range(N), 3) if sum(L) % g == need]
    Es = np.array(list(itertools.combinations_with_replacement(range(N), 3)))
    best = []; allsums = []
    for eps in (0.06, 0.1, 0.15, 0.2, 0.25):
        for h in range(N):
            for L in Ls:
                ne_all = zmin(L[None, :, None] + Es[:, None, :] + h, N)
                se = np.linalg.svd(ce0[None] * eps**ne_all[:, None], compute_uv=False)[..., ::-1]
                okc = np.where(np.all(np.abs(np.median(np.log(se), axis=1) - y_obs) < np.log(3), axis=1))[0]
                if len(okc) == 0: continue
                nn = zmin(L[:, None] + L[None, :] + 2*h, N)
                for ei in okc[:3]:                     # singlet choice barely affects the (hierarchical) charged-lepton rotation
                    t, l, ms = neutrino_stats(nn, ne_all[ei], eps)
                    best.append((t, l, eps, h, L.tolist(), Es[ei].tolist(), float(np.median(ms)) if len(ms) else None)); allsums += list(ms)
    best.sort(key=lambda x: (-x[0], -x[1]))
    print(f"Z_{N}: {len(best)} charged-lepton-compatible assignments tested")
    for b in best[:4]:
        print(f"    tight {b[0]*100:.2f}%  loose {b[1]*100:.1f}%  eps {b[2]} h {b[3]} L {b[4]} e {b[5]}  median sum m_nu (loose hits) {b[6]}")
    if allsums: print(f"    all loose hits: sum m_nu median {np.median(allsums):.3f} eV, 10-90% range {np.percentile(allsums,10):.3f}-{np.percentile(allsums,90):.3f} eV")
