"""Iteration 105: hemisphere test toward/away from (l,b)=(260,+12) (PREREG_105.md)."""
import numpy as np, pandas as pd, os, sys
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar
from astropy.coordinates import SkyCoord; import astropy.units as u
here = os.path.dirname(os.path.abspath(__file__))
ZCUT = float(os.environ.get("ZCUT", 0.05)); NAX = int(os.environ.get("NAX", 300))
PP = "/home/claude/pantheonplussh0es/datarelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/"
sn = pd.read_csv(PP + "Pantheon+SH0ES.dat", sep=r"\s+")
raw = np.loadtxt(PP + "Pantheon+SH0ES_STAT+SYS.cov"); N = int(raw[0]); cov = raw[1:].reshape(N, N)
m = (sn["zHD"] > ZCUT).values; sn = sn[m].reset_index(drop=True); cov = cov[np.ix_(m, m)]
g = SkyCoord(ra=sn.RA.values*u.deg, dec=sn.DEC.values*u.deg).galactic
lr, br = np.radians(g.l.deg), np.radians(g.b.deg)
nh = np.c_[np.cos(br)*np.cos(lr), np.cos(br)*np.sin(lr), np.sin(br)]
z, zh, mb = sn.zHD.values, sn.zHEL.values, sn.m_b_corr.values
zg = np.linspace(0, 2.5, 5001)
def dm_flat(Om): E = np.sqrt(Om*(1+zg)**3 + 1 - Om); return cumulative_trapezoid(1/E, zg, initial=0)
def dm_bare():
    Ob, Ok = 0.049, 0.951; E = np.sqrt(Ob*(1+zg)**3 + Ok*(1+zg)**2); chi = cumulative_trapezoid(1/E, zg, initial=0)
    return np.sinh(np.sqrt(Ok)*chi)/np.sqrt(Ok)
OMS = np.linspace(0.05, 0.8, 151); DMF = {round(o, 4): dm_flat(o) for o in OMS}; DMB = dm_bare()
def c2_of(D, sel, Ci):
    r = mb[sel] - 5*np.log10((1 + zh[sel])*np.interp(z[sel], zg, D)); B = (Ci @ r).sum(); return r @ Ci @ r - B**2/Ci.sum()
def hemi(sel):
    Ci = np.linalg.inv(cov[np.ix_(sel, sel)])
    cs = np.array([c2_of(DMF[round(o, 4)], sel, Ci) for o in OMS]); i = np.argmin(cs)
    ok = OMS[cs <= cs[i] + 1]; om, err = OMS[i], (ok.max() - ok.min())/2
    cb = c2_of(DMB, sel, Ci)
    return dict(n=sel.sum(), om=om, err=max(err, 0.005), pen=(cb - cs[i])/sel.sum(), dpen=cb - cs[i])
def uv(l, b): l, b = np.radians(l), np.radians(b); return np.array([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)])
def stats(ax):
    t = nh @ ax > 0; T, A = hemi(t), hemi(~t)
    return T, A, T["pen"] - A["pen"], T["om"] - A["om"], np.hypot(T["err"], A["err"])
T, A, SA, SB, eB = stats(uv(260, 12))
rng = np.random.default_rng(11); RA_, RB_ = [], []
for _ in range(NAX):
    v = rng.normal(size=3); v /= np.linalg.norm(v); _, _, sa, sb, _ = stats(v); RA_.append(sa); RB_.append(sb)
RA_, RB_ = np.array(RA_), np.array(RB_)
pA = np.mean(np.abs(RA_) >= abs(SA)); pB = np.mean(np.abs(RB_) >= abs(SB))
out = [f"Iteration 105 -- z > {ZCUT}, axis (l,b)=(260,+12), {NAX} random axes",
       f"PART A (no dark energy, no dark matter; ordinary matter only):",
       f"  toward: N={T['n']}, bare-universe penalty vs best LCDM = {T['dpen']:.1f} ({T['pen']:.3f} per SN)",
       f"  away:   N={A['n']}, bare-universe penalty vs best LCDM = {A['dpen']:.1f} ({A['pen']:.3f} per SN)",
       f"  difference per SN {SA:+.3f}; random axes give |diff| this large in {pA:.1%} of cases (spread {RA_.std():.3f})",
       f"PART B (dark energy, LCDM): matter share toward {T['om']:.3f} +/- {T['err']:.3f}, away {A['om']:.3f} +/- {A['err']:.3f}",
       f"  -> dark-energy share toward {1-T['om']:.3f}, away {1-A['om']:.3f}; difference {SB:+.3f} +/- {eB:.3f} ({SB/eB:+.1f} sigma);"
       f" random axes give |diff| this large in {pB:.1%} of cases"]
txt = "\n".join(out); print(txt); open(os.path.join(here, f"iter105_z{ZCUT}.txt"), "w").write(txt + "\n")
