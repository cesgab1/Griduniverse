"""Fence version: sections re-stretched by different amounts, no centre needed.
Prediction: the expansion rate differs by direction on the sky. Test: Pantheon+SH0ES, full covariance,
H0 fitted separately in 12 sky patches (equal-area-ish: 3 declination bands x 4 RA quadrants) and as a dipole."""
import numpy as np, pandas as pd
from scipy.integrate import quad
D = "/home/claude/cobaya_packages/data/sn_data/PantheonPlus/"
d = pd.read_csv(D + "Pantheon+SH0ES.dat", sep=r"\s+"); n = len(d)
C = np.loadtxt(D + "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1).reshape(n, n)
c = 299792.458; Om = 0.315; cal = d.IS_CALIBRATOR.values == 1; z = d.zHD.values
def mu1(zz, zh): return 5*np.log10((1+zh)*c*quad(lambda x: 1/np.sqrt(Om*(1+x)**3+1-Om), 0, zz)[0]) + 25
ra, dec = np.radians(d.RA.values), np.radians(d.DEC.values)
nhat = np.c_[np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)]
for zlo, zhi in ((0.0233, 0.15), (0.15, 2.3)):
    use = cal | ((z > zlo) & (z < zhi)); idx = np.where(use)[0]; Ci = np.linalg.inv(C[np.ix_(idx, idx)])
    y = np.array([d.m_b_corr.values[i] - (d.CEPH_DIST.values[i] if cal[i] else mu1(z[i], d.zHEL.values[i])) for i in idx])
    hub = ~cal[idx]
    # patches
    band = np.digitize(np.sin(dec[idx]), [-1/3, 1/3]); quadr = (np.degrees(ra[idx]) // 90).astype(int) % 4
    patch = band*4 + quadr
    P = 12; X = np.zeros((len(idx), 1+P)); X[:, 0] = 1
    for k in range(len(idx)):
        if hub[k]: X[k, 1+patch[k]] = -5
    keep = [0] + [1+p for p in range(P) if (hub & (patch == p)).sum() >= 5]
    X = X[:, keep]; F = X.T @ Ci @ X; p = np.linalg.solve(F, X.T @ Ci @ y); cov = np.linalg.inv(F)
    H = 10**p[1:]; e = H*np.log(10)*np.sqrt(np.diag(cov)[1:])
    # chi2 of 'all patches equal'
    L = np.log10(H); Cl = cov[1:, 1:]; w = np.linalg.solve(Cl, np.ones(len(L))); mean = (w @ L)/w.sum()
    chi = (L-mean) @ np.linalg.solve(Cl, L-mean)
    print(f"z {zlo}-{zhi}: {hub.sum()} SNe in {len(L)} sky patches; H0 per patch: " + ", ".join(f"{h:.1f}±{s:.1f}" for h, s in zip(H, e)))
    print(f"   spread vs 'same everywhere': chi2 = {chi:.1f} for {len(L)-1} dof;  patch-to-patch rms {np.std(H):.2f}")
    # dipole: log10 H0 = a + b . nhat
    X2 = np.zeros((len(idx), 4)); X2[:, 0] = 1
    X2[hub, 1:] = -5*nhat[idx][hub]
    X3 = np.c_[X2, np.where(hub, -5.0, 0)]
    F = X3.T @ Ci @ X3; q = np.linalg.solve(F, X3.T @ Ci @ y); cq = np.linalg.inv(F)
    b = q[1:4]; amp = np.linalg.norm(b)*np.log(10)*100; eamp = np.sqrt(np.trace(cq[1:4, 1:4])/3)*np.log(10)*100
    print(f"   dipole in expansion rate: {amp:.2f}% (per-axis error {eamp:.2f}%), 95% upper limit about {amp+2*eamp:.1f}%")
