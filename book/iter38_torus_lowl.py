"""
ITERATION 38: the testable half of the book -- is space a flat cube that wraps around, side L ~ 27.5-38 Gpc (rule constant 1;
up to ~60 Gpc for 2)? Test with the REAL Planck 2018 low-l TT likelihood (l = 2-29).
Method: in a wrapped cube only waves that fit (k = 2 pi n / L) exist. Isotropically averaged C_l:
   continuum: C_l = 4 pi  INT dlnk  P(k) Delta_l(k)^2
   torus    : C_l = 4 pi (2 pi)^3 / L^3  SUM_{n != 0}  P(k_n) Delta_l(k_n)^2 / (4 pi k_n^3)
with CAMB's transfer functions Delta_l(k) (Planck 2018 best fit). Ratio torus/continuum multiplies CAMB's D_l.
Validation: L = 150 Gpc must give ratio ~ 1.
Expectations written BEFORE running:
 T1 L in 27.5-38 Gpc suppresses l = 2-3 (missing longest waves), little effect by l ~ 10.
 T2 Planck's low quadrupole makes this a mild IMPROVEMENT or neutral: Delta chi2 between -3 and +1 across the window.
 T3 Not decisive either way: the power spectrum alone cannot confirm a wrap-around (that needs repeated patterns / direction-
    dependent correlations, as in COMPACT's circle searches, which already require L > 27.5 Gpc for a cube).
 Kill: Delta chi2 > +4 across the whole window would disfavour the book's testable half.
"""
import numpy as np, camb
from scipy.signal import fftconvolve
from cobaya.likelihoods.planck_2018_lowl.TT import TT
like = TT({"packages_path": "/home/claude/cobaya_packages"})
pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=60)
pars.set_accuracy(AccuracyBoost=2, lSampleBoost=50)          # every l sampled
res = camb.get_results(pars)
base = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)["total"][:, 0]
td = res.get_cmb_transfer_data("scalar")
q = td.q; Ls = td.L; D = td.delta_p_l_k[0]                    # temperature source: [l, k]
P = lambda k: 2.1e-9*(k/0.05)**(0.9649 - 1)
lmax = 29
def delta_interp(il, k): return np.interp(np.log(k), np.log(q), D[il])
def ratio(Lgpc):
    L = Lgpc*1e3; kmax = 0.03; nmax = int(kmax*L/(2*np.pi)) + 1
    a = np.arange(-nmax, nmax + 1); c1 = np.bincount(a**2, minlength=nmax**2 + 1)
    c2 = np.rint(fftconvolve(c1, c1)); c3 = np.rint(fftconvolve(c2, c1))[:nmax**2 + 1]      # number of lattice vectors with |n|^2 = m
    m = np.arange(1, nmax**2 + 1); cnt = c3[1:]; sel = cnt > 0; m, cnt = m[sel], cnt[sel]
    kn = 2*np.pi*np.sqrt(m)/L
    kf = np.geomspace(1e-6, kmax, 200000)
    out = {}
    for il, l in enumerate(Ls):
        if l < 2 or l > lmax: continue
        cont = np.trapezoid(P(kf)*delta_interp(il, kf)**2/kf, kf)
        tor = (2*np.pi)**3/L**3*np.sum(cnt*P(kn)*delta_interp(il, kn)**2/(4*np.pi*kn**3))
        out[l] = tor/cont
    return out
def chi2(dl): return -2*like.log_likelihood(dl)
c0 = chi2(base)
lines = ["ITERATION 38: wrap-around cube vs the real Planck low-l TT likelihood (expectations committed before running)", "",
         f"plain (infinite space): -2 lnL = {c0:.2f}", "", " L [Gpc]   ratio l=2   l=3    l=5    l=10   l=20    Delta chi2"]
for Lg in (150, 60, 45, 38, 33, 30, 27.5, 23):
    r = ratio(Lg); d = base.copy()
    for l, v in r.items(): d[l] = base[l]*v
    lines.append(f"  {Lg:6.1f}    {r.get(2, np.nan):.3f}   {r.get(3, np.nan):.3f}  {r.get(5, np.nan):.3f}  {r.get(10, np.nan):.3f}  {r.get(20, np.nan):.3f}   {chi2(d) - c0:+.2f}")
txt = "\n".join(lines); print(txt); open("iter38_torus_lowl.txt", "w").write(txt + "\n")
