"""
R3 (truth T3): one acceleration scale. Fitting each SPARC galaxy with the radial-acceleration-relation (RAR) curve
g_obs = nu(g_bar/a0) g_bar gives roughly the same a0 ~ 1.2e-10 m/s^2 everywhere. Does 'ocean = CDM' (Newtonian, NFW halos by
abundance matching with realistic scatter) make galaxies that LOOK like they share one a0?
For each galaxy: best-fit a0 for (i) the observed curve, (ii) the CDM model curve at the same radii (noise-free, so its spread is
purely intrinsic). Also: the global RAR scatter (all points) for data and model.
"""
import numpy as np, glob
from scipy.optimize import brentq, minimize_scalar
rng = np.random.default_rng(4)
G = 4.30e-6; rhoc = 136.0; conv = 3.086e19/1e6          # (km/s)^2/kpc -> m/s^2 : x 1/conv... g[SI] = g / conv
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def shmr(M): M1, N, b, g = 10**11.59, 0.0351, 1.376, 0.608; return 2*N*M/((M/M1)**-b + (M/M1)**g)
def Mh_of(Ms): return np.exp(brentq(lambda lm: np.log(shmr(np.exp(lm))) - np.log(Ms), np.log(1e8), np.log(1e16)))
def slope_ms(Mh): e = 1e-3; return (np.log(shmr(Mh*(1 + e))) - np.log(shmr(Mh)))/np.log(1 + e)
def conc(M): return 10**(0.905 - 0.101*np.log10(M*0.7/1e12))
def nfw_v2(r, M, dc=0.0):
    R200 = (3*M/(4*np.pi*200*rhoc))**(1/3); c = conc(M)*10**dc; x = r*c/R200; m = lambda x: np.log(1 + x) - x/(1 + x)
    return G*M*m(x)/m(c)/r
def fit_a0(gb, go, w):
    f = lambda la: np.sum(w*(np.log10(nu(gb/10**la)*gb) - np.log10(go))**2)
    return minimize_scalar(f, bounds=(-12.5, -8.5), method="bounded").x
A_obs, A_cdm, R_obs, R_cdm, MS = [], [], [], [], []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    ok = (r > 0) & (V > 0); r, V, eV, Vg, Vd, Vb = r[ok], V[ok], eV[ok], Vg[ok], Vd[ok], Vb[ok]
    vb2 = Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2; ok = vb2 > 0
    if ok.sum() < 6: continue
    r, V, eV, vb2, Vd, Vb = r[ok], V[ok], eV[ok], vb2[ok], Vd[ok], Vb[ok]
    Ms = (0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*r[-1]/G
    if Ms < 1e6: continue
    gb = vb2/r/conv; go = V**2/r/conv
    el = 2*np.sqrt(eV**2 + (0.05*V)**2)/V/np.log(10); w = 1/el**2
    MS.append(Ms); A_obs.append(fit_a0(gb, go, w)); R_obs.extend(np.log10(go/(nu(gb/1.2e-10)*gb)))
    Mh = Mh_of(Ms); dh, dc = rng.normal(0, 0.15)/slope_ms(Mh), rng.normal(0, 0.11)
    gm = (vb2 + nfw_v2(r, Mh*10**dh, dc))/r/conv
    A_cdm.append(fit_a0(gb, gm, w)); R_cdm.extend(np.log10(gm/(nu(gb/1.2e-10)*gb)))
A_obs, A_cdm, R_obs, R_cdm = map(np.array, (A_obs, A_cdm, R_obs, R_cdm))
def rob(x): return np.median(x), 0.7413*np.subtract(*np.percentile(x, [75, 25]))
print(f"R3 one acceleration scale: {len(A_obs)} SPARC galaxies, per-galaxy best-fit a0 (fixed M/L 0.5/0.7)")
for lab, A in (("observed", A_obs), ("CDM halos (abundance matched + scatter, noise-free)", A_cdm)):
    m, s = rob(A); print(f"  {lab:52s} median a0 {10**m:.2e} m/s^2   spread {s:.2f} dex (robust)")
print(f"  RAR residuals about the a0 = 1.2e-10 curve, all points: observed mean {np.mean(R_obs):+.3f}, spread {rob(R_obs)[1]:.3f} dex;"
      f" CDM mean {np.mean(R_cdm):+.3f}, spread {rob(R_cdm)[1]:.3f} dex")
# does the CDM a0 drift with galaxy mass? (a real single scale must not)
print("  trend of fitted a0 with galaxy size (log a0 vs rank of stellar mass, terciles):")
MS = np.log10(MS); q = np.percentile(MS, [0, 33.3, 66.7, 100])
for a, b in zip(q[:-1], q[1:]):
    m = (MS >= a) & (MS <= b)
    print(f"    M* 10^{a:.1f}-10^{b:.1f}: observed {10**np.median(A_obs[m]):.2e}   CDM {10**np.median(A_cdm[m]):.2e}")
