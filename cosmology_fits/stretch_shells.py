"""
Is our patch of the grid stretching faster than the universe as a whole?
If regions stretch at different rates, the expansion rate measured nearby (where our patch dominates) differs from farther out:
H0 inferred from supernovae would step down with distance. Test: Pantheon+SH0ES (1701 light curves, full stat+sys covariance),
42 supernovae in galaxies with Cepheid distances fix the brightness scale; H0 fitted separately in distance shells (GLS).
Shape of distance vs redshift beyond ~z 0.1 from flat LCDM with Om = 0.315 (as in SH0ES).
"""
import numpy as np, pandas as pd, json
from scipy.integrate import quad
D = "/home/claude/cobaya_packages/data/sn_data/PantheonPlus/"
d = pd.read_csv(D + "Pantheon+SH0ES.dat", sep=r"\s+")
n = len(d); C = np.loadtxt(D + "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1).reshape(n, n)
c = 299792.458; Om = 0.315
def dl_unit(z, zhel):   # luminosity distance for H0 = 1 (Mpc * km/s/Mpc)
    I = quad(lambda x: 1/np.sqrt(Om*(1+x)**3 + 1 - Om), 0, z)[0]
    return (1 + zhel) * c * I
cal = d.IS_CALIBRATOR.values == 1
edges = [0.01, 0.023, 0.05, 0.1, 0.15, 0.3, 0.6, 2.3]
z = d.zHD.values
use = cal | (z > edges[0])
idx = np.where(use)[0]; Ci = np.linalg.inv(C[np.ix_(idx, idx)])
y = d.m_b_corr.values[idx]; nsh = len(edges) - 1
X = np.zeros((len(idx), 1 + nsh)); X[:, 0] = 1                      # column 0: absolute magnitude M
for k, i in enumerate(idx):
    if cal[i]: y[k] -= d.CEPH_DIST.values[i]
    else:
        y[k] -= 5*np.log10(dl_unit(z[i], d.zHEL.values[i])) + 25     # m - 5log10(dL(H0=1)) - 25 = M - 5 log10 H0
        s = np.searchsorted(edges, z[i]) - 1; X[k, 1 + s] = -5.0    # parameter: log10 H0 in shell s
F = X.T @ Ci @ X; p = np.linalg.solve(F, X.T @ Ci @ y); cov = np.linalg.inv(F)
res = []
print(f"brightness scale from {cal.sum()} Cepheid-calibrated supernovae: M = {p[0]:.3f} +- {np.sqrt(cov[0,0]):.3f}")
for s in range(nsh):
    H = 10**p[1+s]; e = H*np.log(10)*np.sqrt(cov[1+s, 1+s]); nn = int(((z > edges[s]) & (z <= edges[s+1]) & ~cal).sum())
    dist = c*np.sqrt(edges[s]*edges[s+1])/70
    print(f"  shell z {edges[s]:.3f}-{edges[s+1]:.3f} (~{dist:5.0f} Mpc, {nn:4d} SNe): H0 = {H:6.2f} +- {e:.2f}")
    res.append(dict(z0=edges[s], z1=edges[s+1], n=nn, H0=H, err=e))
# step test: shells below vs above z = 0.05 (~210 Mpc) and 0.1
for zc in (0.05, 0.1):
    lo = [s for s in range(nsh) if edges[s+1] <= zc]; hi = [s for s in range(nsh) if edges[s] >= zc]
    def comb(S):
        w = np.zeros(1+nsh); w[[1+s for s in S]] = 1/len(S)
        # weighted by inverse variance sub-block
        sub = np.ix_([1+s for s in S], [1+s for s in S]); Wi = np.linalg.inv(cov[sub]); one = np.ones(len(S))
        wt = Wi @ one / (one @ Wi @ one); v = np.zeros(1+nsh); v[[1+s for s in S]] = wt; return v
    a, b = comb(lo), comb(hi); diff = (a - b) @ p; ediff = np.sqrt((a - b) @ cov @ (a - b))
    print(f"  nearby (z<{zc}) vs farther: local expansion faster by {(10**diff-1)*100:+.2f}% +- {np.log(10)*ediff*100:.2f}%")
json.dump(res, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/stretch_shells.json", "w"), indent=1)
print(f"\nNeeded to explain the tension this way: our patch ~{(73.0/67.9-1)*100:.1f}% faster than the universe as a whole (73.0 vs 67.9)")
