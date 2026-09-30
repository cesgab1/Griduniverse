# Squeezed-layer test: does effective a0 rise with environment (crowded vs isolated, KiDS) or with a galaxy's own depth (SPARC)?
import numpy as np, glob, os
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr
G = 4.52e-30; pcm = 3.086e16; D = "/home/claude/kids/"
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def fit(gb, go, C, gmin=0):
    m = gb > gmin; gb, go, C = gb[m], go[m], C[np.ix_(m, m)]; Ci = np.linalg.inv(C)
    chi = lambda la: (lambda r: r @ Ci @ r)(go - gb*nu(gb/10**la))
    b = minimize_scalar(chi, bounds=(-11.5, -8), method="bounded")
    la = np.linspace(b.x-1, b.x+1, 2001); c2 = np.array([chi(x) for x in la]); ok = la[c2 < b.fun+1]
    return 10**b.x, 10**ok.min(), 10**ok.max(), b.fun, m.sum()
def bins(base, cov):
    c = np.loadtxt(D + cov); keys = np.unique(c[:, 0]); out = []
    for k, key in enumerate(keys):
        d = np.loadtxt(D + base + f"{k+1}.txt"); gb, esd, bias = d[:, 0], d[:, 1], d[:, 4]; n = len(gb)
        blk = c[(c[:, 0] == key) & (c[:, 1] == key)]
        out.append((gb, 4*G*esd/bias*pcm, (blk[:, 4]/blk[:, 6]).reshape(n, n)*(4*G*pcm)**2))
    return out
iso = bins("Fig-9_RAR-KiDS-isolated_Massbin-", "Fig-9_RAR-KiDS-isolated_Massbins_covmatrix.txt")
alls = bins("Fig-A4_RAR-KiDS-all_Massbin-", "Fig-A4_RAR-KiDS-all_Massbins_covmatrix.txt")
labels = ["10^8.5-10.3", "10^10.3-10.6", "10^10.6-10.8", "10^10.8-11.0"]
print("KiDS: a0 for isolated galaxies vs all galaxies (incl. groups/crowded), same stellar-mass bins")
print("Caution: around non-isolated galaxies the neighbours' own mass adds lensing at large radii (weak pull) -> inner-only fit also shown")
for gmin, lab in [(0, "all radii"), (1e-12, "inner only (g_bar > 1e-12)")]:
    print(f"\n  {lab}")
    for L, a, b in zip(labels, iso, alls):
        fi, fa = fit(*a, gmin), fit(*b, gmin)
        print(f"  M* {L}: isolated a0 {fi[0]:.2e} ({fi[1]:.2e}-{fi[2]:.2e})   all {fa[0]:.2e} ({fa[1]:.2e}-{fa[2]:.2e})   ratio all/isolated {fa[0]/fi[0]:.2f}")
# SPARC: per-galaxy a0 vs depth (flat rotation speed)
Yd, Yb, kpc = 0.5, 0.7, 3.0857e19; res = []
for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(f); d = d[None] if d.ndim == 1 else d
    r, V, eV, Vg, Vd, Vb = d.T[:6]
    gb = (Vg*np.abs(Vg) + Yd*Vd**2 + Yb*Vb**2)*1e6/(r*kpc); go = V**2*1e6/(r*kpc)
    ok = (V > 0) & (eV/V < 0.1) & (gb > 0)
    if ok.sum() < 5: continue
    e = 2*eV[ok]/V[ok]/np.log(10)
    chi = lambda la: np.sum(((np.log10(go[ok]) - np.log10(gb[ok]*nu(gb[ok]/10**la)))/np.sqrt(e**2+0.04**2))**2)
    b = minimize_scalar(chi, bounds=(-11.5, -8.5), method="bounded")
    res.append((os.path.basename(f)[:-11], np.median(V[ok][-3:]), b.x))
v = np.array([x[1] for x in res]); la = np.array([x[2] for x in res])
rho, p = spearmanr(v, la); sl = np.polyfit(np.log10(v), la, 1)[0]
print(f"\nSPARC: {len(res)} galaxies, per-galaxy a0 vs outer rotation speed: Spearman ρ = {rho:+.2f} (p = {p:.2g}); log a0 ∝ {sl:+.2f} log v")
for lo, hi in [(0, 80), (80, 150), (150, 220), (220, 400)]:
    m = (v >= lo) & (v < hi)
    print(f"   v {lo}-{hi} km/s: {m.sum():3d} galaxies, median a0 {10**np.median(la[m]):.2e}  (spread {np.std(la[m]):.2f} dex)")
print("Clusters (σ ~ 1000 km/s, v_c ~ 1400): behave as a0 ≈ 2.0e-9 (Tian+2024)")
rng = np.random.default_rng(0); bs = []
for _ in range(2000):
    i = rng.integers(0, len(v), len(v)); bs.append(np.polyfit(np.log10(v[i]), la[i], 1)[0])
print(f"   slope bootstrap: {sl:+.2f} ± {np.std(bs):.2f}; extrapolated to clusters (v 1400): ×{(1400/np.median(v))**sl:.1f};  slope needed for ×18: {np.log10(18)/np.log10(1400/np.median(v)):.2f}")
