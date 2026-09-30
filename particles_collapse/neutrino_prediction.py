"""
Out-of-sample test of the Z_N flavour-charge solutions: neutrinos.
The charges were fitted to quark and charged-lepton masses + CKM only. Neutrino masses come from the Weinberg operator (L_i H)(L_j H):
  exponent n_ij = Z_N-minimal power of (L_i + L_j + 2h);  m_nu ~ c_ij eps^n_ij  (symmetric, O(1) random c_ij)
PMNS = U_e^dagger U_nu. Compared with global-fit values (NuFIT-like): sin^2 th12 = 0.307, sin^2 th23 = 0.47-0.57, sin^2 th13 = 0.022,
Delta m^2_21 / Delta m^2_31 = 0.0296 (normal ordering).
"""
import numpy as np, json, glob
rng = np.random.default_rng(99)
def zmin(n, N): return np.minimum(n % N, N - n % N)
target = dict(s12=0.307, s23=0.52, s13=0.0220, r=0.0296)
for f in sorted(glob.glob("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/flav_ZN_*_w.json")):
    N = int(f.split("_")[-2]); d = json.load(open(f)); ch = np.array(d["charges"]); eps = d["eps"]
    L, e, h = ch[9:12], ch[12:15], ch[15]
    ne = zmin(L[:, None] + e[None, :] + h, N); nn = zmin(L[:, None] + L[None, :] + 2*h, N)
    res = []
    for _ in range(3000):
        ce = np.exp(rng.uniform(np.log(.5), np.log(2), (3, 3))) * np.exp(2j*np.pi*rng.random((3, 3)))
        cn = np.exp(rng.uniform(np.log(.5), np.log(2), (3, 3))) * np.exp(2j*np.pi*rng.random((3, 3))); cn = (cn + cn.T) / 2
        Ue, se, _ = np.linalg.svd(ce * eps**ne)
        Un, sn, _ = np.linalg.svd(cn * eps**nn)          # Takagi ~ SVD for masses and left vectors
        # order: charged leptons e, mu, tau = ascending; neutrinos m1 < m2 < m3
        Ue = Ue[:, ::-1]; Un = Un[:, ::-1]; sn = sn[::-1]
        U = np.abs(Ue.conj().T @ Un)**2
        s13 = U[0, 2]; s12 = U[0, 1] / max(1 - s13, 1e-12); s23 = U[1, 2] / max(1 - s13, 1e-12)
        r = (sn[1]**2 - sn[0]**2) / (sn[2]**2 - sn[0]**2)
        res.append((s12, s23, s13, r))
    res = np.array(res); med = np.median(res, axis=0)
    # fraction of O(1) draws landing within the measured 3-sigma-ish windows
    ok = (res[:, 0] > 0.27) & (res[:, 0] < 0.35) & (res[:, 1] > 0.41) & (res[:, 1] < 0.62) & (res[:, 2] > 0.019) & (res[:, 2] < 0.025) & (res[:, 3] > 0.02) & (res[:, 3] < 0.045)
    ok_ang = (res[:, 0] > 0.2) & (res[:, 0] < 0.45) & (res[:, 1] > 0.3) & (res[:, 1] < 0.7) & (res[:, 2] > 0.005) & (res[:, 2] < 0.08)
    print(f"Z_{N:<2d} (lepton doublet charges {L}, eps = {eps:.3f}, Weinberg exponents {nn.tolist()}):")
    print(f"     median sin^2 th12 = {med[0]:.3f} (meas 0.307), sin^2 th23 = {med[1]:.3f} (0.47-0.57), sin^2 th13 = {med[2]:.4f} (0.022), dm21/dm31 = {med[3]:.3f} (0.030)")
    print(f"     share of O(1) draws matching all four within ~3 sigma: {ok.mean()*100:.1f}%;  matching the three angles loosely (factor ~2): {ok_ang.mean()*100:.1f}%")
