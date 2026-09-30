"""
C. The missing ingredient: three right-handed neutrinos N_R with Z_N flavour charges R_j (grid-allowed: they carry no SM charge,
   so they add nothing to the SU(3)^2 and SU(2)^2 anomaly rule).
   Dirac coupling  Y_D ij = c eps^min(L_i + R_j mod N),   Majorana mass  M_R ij = M_grid c eps^min(R_i + R_j mod N),   M_grid = 1/l_d
   Light neutrinos: m_nu = (v^2/2) Y_D M_R^-1 Y_D^T  (type-I seesaw).  Everything else as in the lepton census.
   For each (L, R) we draw 200 sets of O(1) coefficients and count how often ALL of these come out right:
     loose: sin^2 th12 [.2,.45], th23 [.35,.65], th13 [.01,.05], dm21/dm31 [.01,.08], heaviest mass m3 in [0.035, 0.1] eV
     tight: [.27,.35], [.41,.62], [.019,.025], [.024,.036], m3 in [0.045, 0.07] eV
   Baseline (Weinberg operator, no N_R): best case 0.2% loose, 0% tight.
"""
import numpy as np, itertools, json, glob, sys
from math import gcd
rng = np.random.default_rng(5)
N = int(sys.argv[1]); v = 246.22; Mg = 1.220890e19 / 0.603
ND0, ND = 24, 200
y_obs = np.log(np.array([0.4866e-3, 0.1027, 1.746]) * np.sqrt(2) / v)
def cO1(n, sym=False):
    c = np.exp(rng.uniform(np.log(.5), np.log(2), (n, 3, 3))) * np.exp(2j*np.pi*rng.random((n, 3, 3)))
    return (c + np.transpose(c, (0, 2, 1))) / 2 if sym else c
ce0, ce, cD, cM = cO1(ND0), cO1(ND), cO1(ND), cO1(ND, True)
def zmin(n): return np.minimum(n % N, N - n % N)
Q = np.array(json.load(open(glob.glob(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/flav_ZN_{N}_nu.json")[0]))["charges"][0:3])
g = gcd(9, N); need = (-3 * Q.sum()) % g
Ls = [np.array(L) for L in itertools.combinations_with_replacement(range(N), 3) if sum(L) % g == need]
Rs = np.array(list(itertools.combinations_with_replacement(range(N), 3)))
Es = np.array(list(itertools.combinations_with_replacement(range(N), 3)))
results = []
for eps in (0.15, 0.2, 0.25):
    for L in Ls:
        # charged leptons: best singlet assignment for this L (h = 0; a Higgs-charge shift is absorbed into e and R)
        ne_all = zmin(L[None, :, None] + Es[:, None, :])
        se = np.linalg.svd(ce0[None] * eps**ne_all[:, None], compute_uv=False)[..., ::-1]
        dev = np.max(np.abs(np.median(np.log(se), axis=1) - y_obs), axis=1)
        if dev.min() > np.log(3): continue
        ne = ne_all[np.argmin(dev)]
        Ue, _, _ = np.linalg.svd(ce * eps**ne); Ue = Ue[:, :, ::-1]
        nD = zmin(L[None, :, None] + Rs[:, None, :])                    # (nR, 3, 3)
        nM = zmin(Rs[:, :, None] + Rs[:, None, :])
        YD = cD[None] * eps**nD[:, None]; MR = Mg * cM[None] * eps**nM[:, None]
        mnu = (v**2 / 2) * YD @ np.linalg.inv(MR) @ np.transpose(YD, (0, 1, 3, 2)) * 1e9     # eV
        Un, sn, _ = np.linalg.svd(mnu); Un = Un[..., ::-1]; sn = sn[..., ::-1]
        P = np.abs(np.einsum("dji,rdjl->rdil", Ue.conj(), Un))**2
        s13 = P[..., 0, 2]; s12 = P[..., 0, 1] / (1 - s13 + 1e-12); s23 = P[..., 1, 2] / (1 - s13 + 1e-12)
        r = (sn[..., 1]**2 - sn[..., 0]**2) / (sn[..., 2]**2 - sn[..., 0]**2 + 1e-300); m3 = sn[..., 2]
        loose = (s12 > .2) & (s12 < .45) & (s23 > .35) & (s23 < .65) & (s13 > .01) & (s13 < .05) & (r > .01) & (r < .08) & (m3 > .035) & (m3 < .1)
        tight = (s12 > .27) & (s12 < .35) & (s23 > .41) & (s23 < .62) & (s13 > .019) & (s13 < .025) & (r > .024) & (r < .036) & (m3 > .045) & (m3 < .07)
        pl, pt = loose.mean(1), tight.mean(1)
        for k in np.argsort(-pl)[:3]:
            if pl[k] == 0: continue
            hit = loose[k]
            m1 = sn[k, hit, 0]; m2 = sn[k, hit, 1]; m3h = sn[k, hit, 2]
            Up = np.einsum("dji,djl->dil", Ue[hit].conj(), Un[k, hit])
            mbb = np.abs(np.einsum("di,di->d", Up[:, 0, :]**2, sn[k, hit]))
            MRdiag = np.median(np.sort(np.abs(np.linalg.eigvals(MR[k])), axis=1), axis=0)
            results.append(dict(eps=eps, L=L.tolist(), R=Rs[k].tolist(), p_loose=float(pl[k]), p_tight=float(pt[k]),
                                sum_mnu=float(np.median(m1 + m2 + m3h)), m1=float(np.median(m1)), mbb=float(np.median(mbb)),
                                MR_GeV=[float(x) for x in MRdiag]))
results.sort(key=lambda x: (-x["p_tight"], -x["p_loose"]))
print(f"Z_{N}: {len(Ls)} lepton-doublet choices x {len(Rs)} right-handed-neutrino choices x 3 eps")
for r_ in results[:6]:
    print(f"   p(tight) {r_['p_tight']*100:5.1f}%  p(loose) {r_['p_loose']*100:5.1f}%  eps {r_['eps']}  L {r_['L']}  R {r_['R']}  | sum m_nu {r_['sum_mnu']:.3f} eV,"
          f" m_lightest {r_['m1']*1e3:.1f} meV, m_bb {r_['mbb']*1e3:.1f} meV, heavy N masses {['%.0e' % x for x in r_['MR_GeV']]} GeV")
good = [r_ for r_ in results if r_["p_loose"] >= 0.05]
if good:
    S = np.array([r_["sum_mnu"] for r_ in good]); B = np.array([r_["mbb"] for r_ in good]); M1 = np.array([r_["m1"] for r_ in good])
    print(f"   {len(good)} assignments reach >= 5% (loose). Their predictions: sum m_nu {np.percentile(S,10):.3f}-{np.percentile(S,90):.3f} eV (median {np.median(S):.3f}),"
          f" m_bb {np.percentile(B,10)*1e3:.1f}-{np.percentile(B,90)*1e3:.1f} meV, lightest {np.percentile(M1,10)*1e3:.1f}-{np.percentile(M1,90)*1e3:.1f} meV")
json.dump(results[:500], open(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/seesaw_Z{N}.json", "w"))
