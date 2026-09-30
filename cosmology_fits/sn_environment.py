"""
Is the supernova 'tape measure' bent by how taut the grid is where each supernova explodes?
Tension picture: local G ~ tautness. White-dwarf supernova energy tracks the Chandrasekhar mass ~ G^(-3/2):
   stronger local G -> dimmer:  delta m = 1.63 * delta G / G.
Two proxies for tautness at the explosion site (Pantheon+SH0ES, standard corrections already applied, incl. the usual host-mass step):
   P1 = host galaxy stellar mass (log10 M*)
   P2 = local pull ~ M* / r^2 at the projected distance of the supernova from its host's centre
Tests: (1) do calibrator hosts (Cepheid galaxies) differ from distant 'Hubble-flow' hosts?
       (2) fit H0 jointly with a tautness-dependent G:  m = M + mu(z, H0) + 1.63 k (P - P_ref).  Does H0 move toward 68?
"""
import numpy as np, pandas as pd
from scipy.integrate import quad
D = "/home/claude/cobaya_packages/data/sn_data/PantheonPlus/"
d = pd.read_csv(D + "Pantheon+SH0ES.dat", sep=r"\s+"); n = len(d)
C = np.loadtxt(D + "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1).reshape(n, n)
c = 299792.458; Om = 0.315; z = d.zHD.values; cal = d.IS_CALIBRATOR.values == 1
hf = (~cal) & (z > 0.0233) & (z < 0.15)
def dl1(zz, zh): return (1+zh)*c*quad(lambda x: 1/np.sqrt(Om*(1+x)**3+1-Om), 0, zz)[0]
dist = np.array([10**((d.CEPH_DIST.values[i]-25)/5) if cal[i] else dl1(z[i], d.zHEL.values[i])/73/(1+d.zHEL.values[i])**2 for i in range(n)])  # Mpc (angular-diameter for HF)
lm = d.HOST_LOGMASS.values
rproj = np.radians(d.HOST_ANGSEP.values/3600) * dist * 1e3  # kpc
good_r = (d.HOST_ANGSEP.values > 0) & (rproj > 0.3)
pull = lm - 2*np.log10(np.clip(rproj, 0.3, None))
print("1. Do the two groups of galaxies differ?")
for lab, P, ok in (("host stellar mass log M*", lm, lm > 0), ("local pull log(M*/r^2)", pull, (lm > 0) & good_r)):
    a, b = P[cal & ok], P[hf & ok]
    print(f"   {lab}: calibrator hosts median {np.median(a):.2f} (N={len(a)}), distant hosts median {np.median(b):.2f} (N={len(b)}); difference {np.median(a)-np.median(b):+.2f}")
idx = np.where(cal | hf)[0]; Ci = np.linalg.inv(C[np.ix_(idx, idx)])
y = np.array([d.m_b_corr.values[i] - (d.CEPH_DIST.values[i] if cal[i] else 5*np.log10(dl1(z[i], d.zHEL.values[i]))+25) for i in idx])
isf = hf[idx].astype(float)
def fit(P=None, ok=None):
    cols = [np.ones(len(idx)), -5*isf]
    if P is not None:
        Pc = np.where(ok[idx], P[idx] - np.median(P[idx][ok[idx]]), 0.0); cols.append(1.63*Pc)
    X = np.c_[tuple(cols)]; F = X.T @ Ci @ X; p = np.linalg.solve(F, X.T @ Ci @ y); cv = np.linalg.inv(F)
    r = y - X @ p; return p, cv, r @ Ci @ r
p0, c0, chi0 = fit()
H = 10**p0[1]; print(f"\n2. Standard fit: H0 = {H:.2f} +- {H*np.log(10)*np.sqrt(c0[1,1]):.2f}   ({len(idx)} light curves)")
for lab, P, ok in (("host stellar mass", lm, lm > 0), ("local pull M*/r^2", pull, (lm > 0) & good_r)):
    p, cv, chi = fit(P, ok); H = 10**p[1]
    print(f"   with G tied to {lab}: k = dG/G per dex = {p[2]*100:+.2f}% +- {np.sqrt(cv[2,2])*100:.2f}%  -> H0 = {H:.2f} +- {H*np.log(10)*np.sqrt(cv[1,1]):.2f}; fit improves by {chi0-chi:.1f}")
# what would be needed
dP = np.median(lm[cal & (lm > 0)]) - np.median(lm[hf & (lm > 0)])
need_mag = 5*np.log10(73.0/67.9)
print(f"\n3. To bring H0 to 67.9 the calibrator supernovae must be effectively {need_mag:.2f} mag brighter than assumed relative to distant ones:")
print(f"   = local G weaker by {need_mag/1.63*100:.0f}% in calibrator hosts; with a mass difference of {dP:+.2f} dex that needs k = {-need_mag/1.63/dP*100 if dP else float('nan'):+.0f}% per dex")
