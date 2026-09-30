"""
What neutrino mixing needs: compare, with the same O(1)-coefficient statistics, how often each structure reproduces the measured angles.
  1. Anarchy (no flavour structure; best abelian Z_N case)              m_nu = random O(1) symmetric matrix
  2. Non-abelian discrete symmetry (A4/S4 type) -> tri-bimaximal (TBM) mixing at leading order, corrected at order eps:
        m_nu = U_TBM* diag(m1, m2, m3) U_TBM^dagger + eps m3 C,   U_e = exp(i eps H)   (C, H random O(1))
     with the same eps family as the quark sector (0.05-0.2). Masses: normal ordering, m1 ~ 0.
Windows: tight ~3 sigma: sin^2 th12 [.27,.35], th23 [.41,.62], th13 [.019,.025];  loose: [.2,.45], [.35,.65], [.01,.05].
"""
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(3)
ND = 20000
def angles(Ue, Un):
    P = np.abs(np.einsum("dji,djl->dil", Ue.conj(), Un))**2
    s13 = P[:, 0, 2]; return P[:, 0, 1] / (1 - s13), P[:, 1, 2] / (1 - s13), s13
def windows(s12, s23, s13):
    t = (s12 > .27) & (s12 < .35) & (s23 > .41) & (s23 < .62) & (s13 > .019) & (s13 < .025)
    l = (s12 > .2) & (s12 < .45) & (s23 > .35) & (s23 < .65) & (s13 > .01) & (s13 < .05)
    return t.mean(), l.mean()
def cO1(sym=False, herm=False):
    c = np.exp(rng.uniform(np.log(.5), np.log(2), (ND, 3, 3))) * np.exp(2j*np.pi*rng.random((ND, 3, 3)))
    if sym: c = (c + np.transpose(c, (0, 2, 1))) / 2
    if herm: c = (c + np.transpose(c.conj(), (0, 2, 1))) / 2
    return c
def order(M):
    U, s, _ = np.linalg.svd(M); return U[:, :, ::-1], s[:, ::-1]
I3 = np.broadcast_to(np.eye(3), (ND, 3, 3))
Un, sn = order(cO1(sym=True))
t, l = windows(*angles(I3, Un))
print(f"1. anarchy (random O(1) mass matrix):                         tight {t*100:5.2f}%   loose {l*100:5.1f}%")
TBM = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0], [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)], [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)]])
m = np.array([0.0, np.sqrt(7.42e-5), np.sqrt(2.51e-3)])
M0 = TBM.conj() @ np.diag(m) @ TBM.conj().T
print(f"   (exact tri-bimaximal: sin^2 th12 = 1/3, th23 = 1/2, th13 = 0 -> th13 is the one that must come from corrections)")
out = {}
for eps in (0.05, 0.1, 0.15, 0.2):
    C = cO1(sym=True); H = cO1(herm=True)
    Mn = M0[None] + eps * m[2] * C
    Un, sn = order(Mn)
    Ue = np.array([expm(1j * eps * h) for h in H[:4000]]); Unn = Un[:4000]
    s12, s23, s13 = angles(Ue, Unn); t, l = windows(s12, s23, s13)
    r = (sn[:4000, 1]**2 - sn[:4000, 0]**2) / (sn[:4000, 2]**2 - sn[:4000, 0]**2)
    ok_r = np.mean((r > .024) & (r < .036))
    print(f"2. tri-bimaximal + corrections of size eps = {eps:4.2f}:            tight {t*100:5.2f}%   loose {l*100:5.1f}%   "
          f"| median sin^2 th13 {np.median(s13):.4f}, th23 spread (10-90%) {np.percentile(s23,10):.2f}-{np.percentile(s23,90):.2f}, "
          f"splitting ratio kept in range {ok_r*100:.0f}%")
    out[eps] = (t, l)
best = max(out, key=lambda e: out[e][0])
print(f"\n=> best: eps = {best}: {out[best][0]*100:.1f}% tight vs anarchy {windows(*angles(I3, order(cO1(sym=True))[0]))[0]*100:.2f}%")
