# Does a0 grow with lookback time? Fit a0 to KiDS/GAMA weak-lensing RAR (lenses at z ~ 0.1-0.5) and compare to SPARC (z ≈ 0)
import numpy as np
from scipy.optimize import minimize_scalar
G = 4.52e-30; pcm = 3.086e16; D = "/home/claude/kids/"
def load(fn, cov):
    d = np.loadtxt(D + fn); gb, esd, bias = d[:, 0], d[:, 1], d[:, 4]
    c = np.loadtxt(D + cov); n = len(gb)
    return gb, 4*G*esd/bias*pcm, (c[:, 4]/c[:, 6]).reshape(n, n)*(4*G*pcm)**2
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
res={}
for lab, f, cf in [("KiDS isolated (stars+cold gas)","Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt","Fig-4-5-C1_RAR-KiDS-isolated_covmatrix.txt"),
                   ("KiDS isolated + hot gas","Fig-4_RAR-KiDS-isolated_hotgas_Nobins.txt","Fig-4_RAR-KiDS-isolated_hotgas_covmatrix.txt"),
                   ("GAMA isolated","Fig-4-C1_RAR-GAMA-isolated_Nobins.txt","Fig-4-C1_RAR-GAMA-isolated_covmatrix.txt")]:
    gb, go, C = load(f, cf); Ci = np.linalg.inv(C)
    chi = lambda la: (lambda r: r @ Ci @ r)(go - gb*nu(gb/10**la))
    b = minimize_scalar(chi, bounds=(-11, -8.5), method="bounded")
    la = np.linspace(b.x-0.6, b.x+0.6, 1201); c2 = np.array([chi(x) for x in la]); ok = la[c2 < b.fun+1]
    a0, lo, hi = 10**b.x, 10**ok.min(), 10**ok.max()
    res[lab]=(a0,lo,hi)
    print(f"{lab:32s} a0 = {a0:.2e}  (1σ {lo:.2e}–{hi:.2e})  chi2 {b.fun:.1f}/{len(gb)-1}")
sparc=1.15e-10
print(f"\nSPARC (z≈0): {sparc:.2e}")
for z in (0.2, 0.3):
    m=(1.0+1.59*z)/1.0
    print(f"MUSE-DARK growth predicts at lens z={z}: ×{m:.2f} -> {sparc*m:.2e};  a0 ∝ H: ×{np.sqrt(0.31*(1+z)**3+0.69):.2f};  tension-energy rule: ×1.00")
a0,lo,hi=res["KiDS isolated (stars+cold gas)"]
print(f"\nKiDS/SPARC = {a0/sparc:.2f} (+{hi/sparc-a0/sparc:.2f} / -{a0/sparc-lo/sparc:.2f})")

print("\n=== split by galaxy type (discs carry little hot gas, so a0 is cleaner) ===")
for base, cov, names in [("Fig-8_RAR-KiDS-isolated_Colorbin_", "Fig-8_RAR-KiDS-isolated_Colorbins_covmatrix.txt", ("blue (discs)","red (ellipticals)")),
                         ("Fig-8_RAR-KiDS-isolated_Sersicbin_", "Fig-8_RAR-KiDS-isolated_Sersicbins_covmatrix.txt", ("low Sérsic (discs)","high Sérsic (ellipticals)"))]:
    c = np.loadtxt(D + cov); keys = np.unique(c[:, 0])
    for k, key in enumerate(keys):
        d = np.loadtxt(D + base + f"{k+1}.txt"); gb, esd, bias = d[:, 0], d[:, 1], d[:, 4]; n = len(gb)
        go = 4*G*esd/bias*pcm
        blk = c[(c[:, 0] == key) & (c[:, 1] == key)]
        C = (blk[:, 4]/blk[:, 6]).reshape(n, n)*(4*G*pcm)**2; Ci = np.linalg.inv(C)
        chi = lambda la: (lambda r: r @ Ci @ r)(go - gb*nu(gb/10**la))
        b = minimize_scalar(chi, bounds=(-11, -8.5), method="bounded")
        la = np.linspace(b.x-0.8, b.x+0.8, 1601); c2 = np.array([chi(x) for x in la]); ok = la[c2 < b.fun+1]
        print(f"{names[k]:28s} a0 = {10**b.x:.2e} (1σ {10**ok.min():.2e}–{10**ok.max():.2e})  ratio to SPARC {10**b.x/sparc:.2f}  chi2 {b.fun:.1f}/{n-1}")
