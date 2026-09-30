"""
Census of ALL grid-compatible lepton charge assignments for a discrete Z_N flavour symmetry (N = 9, 12, 15, 18).
Fixed from each N's quark fit: quark doublet charges Q (they enter the anomaly rule).
Free: Higgs charge h (a shift of h with compensating u, d shifts leaves every quark exponent unchanged) and eps in {0.06..0.25}.
Anomaly rule: sum(3Q + L) = 0 mod N; a common shift of the quark charges leaves the quark sector unchanged, so only
  sum(L) = -3 sum(Q)  (mod gcd(9, N)) is a real constraint.
Exhaustive over lepton doublets L_i and singlets e_j in Z_N (sorted; generation relabelling is not physical).
Charged leptons: y_ij ~ c eps^min(L_i+e_j+h mod N);  neutrinos (Weinberg): m_ij ~ c eps^min(L_i+L_j+2h mod N), c = O(1) random.
Solution: median e, mu, tau Yukawas within x3 AND sin^2 th12 in [.2,.45], sin^2 th23 in [.35,.65], sin^2 th13 in [.01,.05], dm21/dm31 in [.01,.08].
"""
import numpy as np, json, itertools, glob
from math import gcd
rng = np.random.default_rng(2024)
ND = 24; v = 246.22
y_obs = np.log(np.array([0.4866e-3, 0.1027, 1.746]) * np.sqrt(2) / v)
ce = np.exp(rng.uniform(np.log(.5), np.log(2), (ND, 3, 3))) * np.exp(2j*np.pi*rng.random((ND, 3, 3)))
cn = np.exp(rng.uniform(np.log(.5), np.log(2), (ND, 3, 3))) * np.exp(2j*np.pi*rng.random((ND, 3, 3))); cn = (cn + np.transpose(cn, (0, 2, 1))) / 2
dm31, dm21 = 2.51e-3, 7.42e-5
def zmin(n, N): return np.minimum(n % N, N - n % N)

def scan(N, Q):
    g = gcd(9, N); need = (-3 * Q.sum()) % g
    Ls = [np.array(L) for L in itertools.combinations_with_replacement(range(N), 3) if sum(L) % g == need]
    Es = np.array(list(itertools.combinations_with_replacement(range(N), 3)))
    sols, near, n_charged = [], [], 0
    for eps in (0.06, 0.1, 0.15, 0.2, 0.25):
        for h in range(N):
            for L in Ls:
                nn = zmin(L[:, None] + L[None, :] + 2*h, N)
                Un, sn, _ = np.linalg.svd(cn * eps**nn); Un = Un[:, :, ::-1]; sn = sn[:, ::-1]
                ne = zmin(L[None, :, None] + Es[:, None, :] + h, N)
                Ue, se, _ = np.linalg.svd(ce[None] * eps**ne[:, None]); Ue = Ue[..., ::-1]; se = se[..., ::-1]
                okc = np.all(np.abs(np.median(np.log(se), axis=1) - y_obs) < np.log(3), axis=1)
                n_charged += int(okc.sum())
                if not okc.any(): continue
                idx = np.where(okc)[0]
                Up = np.einsum("edji,djl->edil", Ue[idx].conj(), Un); P = np.abs(Up)**2
                s13 = P[..., 0, 2]; s12 = P[..., 0, 1] / (1 - s13 + 1e-12); s23 = P[..., 1, 2] / (1 - s13 + 1e-12)
                med = np.c_[np.median(s12, 1), np.median(s23, 1), np.median(s13, 1)]
                rm = float(np.median((sn[:, 1]**2 - sn[:, 0]**2) / (sn[:, 2]**2 - sn[:, 0]**2 + 1e-300)))
                ratio = float(np.median(sn[:, 0] / sn[:, 2]))
                dist = (np.abs(np.log(med[:, 0] / .307)) + np.abs(np.log(med[:, 1] / .52)) + np.abs(np.log(med[:, 2] / .022))
                        + abs(np.log(max(rm, 1e-9) / .0296)))
                k0 = int(np.argmin(dist))
                near.append((round(float(dist[k0]), 2), eps, h, L.tolist(), Es[idx[k0]].tolist(), med[k0].round(3).tolist(), round(rm, 3), round(ratio, 3)))
                ok = (med[:, 0] > .2) & (med[:, 0] < .45) & (med[:, 1] > .35) & (med[:, 1] < .65) & (med[:, 2] > .01) & (med[:, 2] < .05) & (0.01 < rm < 0.08)
                if ok.any():
                    m3 = np.sqrt(dm31 / (1 - ratio**2)); m1 = ratio * m3; m2 = np.sqrt(m1**2 + dm21)
                    mbb = np.median(np.abs(np.einsum("edi,di->ed", Up[ok][:, :, 0, :]**2, sn / sn[:, 2:3])), axis=1) * m3
                    for k, ei in enumerate(idx[ok]):
                        sols.append(dict(eps=eps, h=h, L=L.tolist(), e=Es[ei].tolist(), s12=float(med[ok][k, 0]), s23=float(med[ok][k, 1]),
                                         s13=float(med[ok][k, 2]), r=rm, m1_m3=ratio, sum_mnu=float(m1 + m2 + m3), mbb=float(mbb[k])))
    return Ls, Es, sols, near, n_charged

summary = {}
for f in sorted(glob.glob("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/flav_ZN_*_nu.json"), key=lambda x: int(x.split("_")[-2])):
    N = int(f.split("_")[-2]); Q = np.array(json.load(open(f))["charges"][0:3])
    Ls, Es, sols, near, n_charged = scan(N, Q)
    print(f"Z_{N:<2d}: {len(Ls)} doublet x {len(Es)} singlet choices x {N} Higgs charges x 5 eps; charged leptons fit in {n_charged}; neutrinos also fit in {len(sols)}")
    if sols:
        S = np.array([s["sum_mnu"] for s in sols]); mb = np.array([s["mbb"] for s in sols])
        print(f"      predicted sum m_nu: min {S.min():.3f}  median {np.median(S):.3f}  max {S.max():.3f} eV;  below 0.10 eV: {np.mean(S < 0.10)*100:.0f}%;"
              f"  m_bb {mb.min()*1e3:.0f}-{mb.max()*1e3:.0f} meV;  m1/m3 {min(s['m1_m3'] for s in sols):.2f}-{max(s['m1_m3'] for s in sols):.2f}")
    near.sort(key=lambda x: x[0])
    print("      closest (distance, eps, h, L, e, [s12, s23, s13], dm-ratio, m1/m3):"); [print("       ", x) for x in near[:3]]
    summary[N] = dict(n_solutions=len(sols), n_charged=n_charged, solutions=sols[:300], near=near[:20])
json.dump(summary, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/lepton_census.json", "w"))
