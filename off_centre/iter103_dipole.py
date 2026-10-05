"""Iteration 103: direction (dipole) test of the supernova Hubble diagram (PREREG_103.md)."""
import numpy as np, pandas as pd, os
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar
from scipy.stats import chi2 as chi2d
from astropy.coordinates import SkyCoord
import astropy.units as u
here = os.path.dirname(os.path.abspath(__file__))
PP = "/home/claude/pantheonplussh0es/datarelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/"
sn = pd.read_csv(PP + "Pantheon+SH0ES.dat", sep=r"\s+")
raw = np.loadtxt(PP + "Pantheon+SH0ES_STAT+SYS.cov"); N = int(raw[0]); cov = raw[1:].reshape(N, N)
m = (sn["zHD"] > 0.01).values; sn = sn[m].reset_index(drop=True); cov = cov[np.ix_(m, m)]
ra, dec = np.radians(sn.RA.values), np.radians(sn.DEC.values)
nhat = np.c_[np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)]
zg = np.linspace(0, 2.5, 5001)
def DM(Om):
    E = np.sqrt(Om*(1 + zg)**3 + 1 - Om); return cumulative_trapezoid(1/E, zg, initial=0)
def gls(X, y, C):
    Ci = np.linalg.inv(C); F = X.T @ Ci @ X; p = np.linalg.solve(F, X.T @ Ci @ y)
    r = y - X @ p; return p, np.linalg.inv(F), r @ Ci @ r
def vec2lb(v):
    c = SkyCoord(x=v[0], y=v[1], z=v[2], representation_type="cartesian", frame="icrs"); g = c.galactic
    return g.l.deg, g.b.deg
ref = {"CMB": SkyCoord(l=264.0*u.deg, b=48.3*u.deg, frame="galactic"), "quasars": SkyCoord(l=238.2*u.deg, b=28.8*u.deg, frame="galactic")}
out = []
for zcol in ("zCMB", "zHD"):
    z = sn[zcol].values; zhel = sn.zHEL.values; mb = sn.m_b_corr.values
    def mu(Om): return 5*np.log10((1 + zhel)*np.interp(z, zg, DM(Om))) + 25
    Ci = np.linalg.inv(cov)
    def c2(Om):
        r = mb - mu(Om); B = (Ci @ r).sum(); return r @ Ci @ r - B**2/Ci.sum()
    Om = minimize_scalar(c2, bounds=(0.1, 0.6), method="bounded").x
    res = mb - mu(Om)
    out.append(f"== redshift {zcol}: baseline Om {Om:.3f} ==")
    for lab, sel in (("z 0.01-0.05", (z < 0.05)), ("z 0.05-0.15", (z >= 0.05) & (z < 0.15)), ("z > 0.15", z >= 0.15), ("all", z > 0)):
        C = cov[np.ix_(sel, sel)]; y = res[sel]; n = nhat[sel]
        X0 = np.ones((sel.sum(), 1)); _, _, c0 = gls(X0, y, C)
        X1 = np.c_[np.ones(sel.sum()), -5/np.log(10)*n]          # distance x (1 + D.n): mu shift = 5/ln10 * D.n ... sign: faster expansion -> smaller distance
        p, cv, c1 = gls(X1, y, C)
        D = p[1:]; A = np.linalg.norm(D); eA = np.sqrt(D @ cv[1:, 1:] @ D)/A
        dchi = c0 - c1; pval = chi2d.sf(dchi, 3)
        l, b = vec2lb(D); dirc = SkyCoord(l=l*u.deg, b=b*u.deg, frame="galactic")
        sep = {k: dirc.separation(v).deg for k, v in ref.items()}
        # 95% upper limit on amplitude from Monte Carlo of the fitted covariance
        sims = np.random.default_rng(0).multivariate_normal(np.zeros(3), cv[1:, 1:], 4000)
        lim = np.percentile(np.linalg.norm(sims, axis=1), 95)
        out.append(f"  {lab:12s} N={sel.sum():4d}  A = {A*100:5.2f}% +/- {eA*100:.2f}%  toward (l,b)=({l:5.1f},{b:+5.1f})"
                   f"  delta chi2 {dchi:5.2f} (p={pval:.3f})  | angle to CMB dipole {sep['CMB']:5.1f}, to quasar dipole {sep['quasars']:5.1f}"
                   f"  | noise-only 95% amplitude {lim*100:.2f}%")
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter103_dipole.txt"), "w").write(txt + "\n")
