"""
Step 2: can 'cells grow where the pull is weak' produce the galaxy law by itself?
Grid mechanics (string network): effective tension = link tension x (total link length per volume). In 3D, cells of size s have
link length per volume ~ 1/s^2, so G_eff/G = (s/s0)^2  [checked numerically below on random grids of different density].
Cracking rules for how cell size follows the local pull g (strain ~ g):
   Griffith cracking (stored energy vs crack surface):          s ~ g^-2
   crack spacing inversely proportional to strain (thin films): s ~ g^-1
   generic power law s ~ g^-p ; boost = (s/s0)^2 = 1 + (a_c/g)^(2p)
Test on 164 SPARC galaxies (2,803 points): fit a_c for each p; the observed law needs boost ~ sqrt(a0/g) at weak pull, i.e. p = 1/4.
"""
import numpy as np, glob
from scipy.spatial import Voronoi
from scipy.optimize import minimize_scalar, minimize
rng = np.random.default_rng(0)
print("check: link length per volume vs cell size on random 3D grids")
for n in (500, 2000, 8000):
    L = 1.0; P = rng.uniform(0, L, (n, 3)); vor = Voronoi(P); I, J = vor.ridge_points.T
    inside = (np.abs(P[I] - 0.5).max(1) < 0.3) & (np.abs(P[J] - 0.5).max(1) < 0.3)
    dens = np.linalg.norm(P[I] - P[J], axis=1)[inside].sum()/0.6**3; s = n**(-1/3)
    print(f"   cell size {s:.3f}: link length per volume {dens:8.1f}, times s^2 = {dens*s**2:.2f} (constant -> tension ~ 1/s^2, G_eff ~ s^2)")
Yd, Yb, kpc = 0.5, 0.7, 3.0857e19; go, gb, el = [], [], []
for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(f, comments="#"); d = d[None] if d.ndim == 1 else d
    r, V, eV, Vg, Vd, Vb = d.T[:6]
    g1 = (Vg*np.abs(Vg) + Yd*Vd**2 + Yb*Vb**2)*1e6/(r*kpc); g0 = V**2*1e6/(r*kpc)
    ok = (V > 0) & (eV/V < 0.1) & (g1 > 0)
    go += list(g0[ok]); gb += list(g1[ok]); el += list(2*eV[ok]/V[ok]/np.log(10))
go, gb, el = np.log10(go), np.array(gb), np.array(el)
def chi(p, lac): pred = np.log10(gb*(1 + (10**lac/gb)**(2*p))); return np.sum(((go - pred)/np.hypot(el, 0.1))**2)
print("\nfit to 164 galaxies, boost = (cell size / small-cell size)^2 = 1 + (a_c/g)^(2p):")
for p, lab in ((2, "Griffith cracking, s ~ g^-2"), (1, "crack spacing ~ 1/strain, s ~ g^-1"), (0.5, "s ~ g^-1/2"), (0.25, "s ~ g^-1/4 (what the galaxies need)")):
    b = minimize_scalar(lambda x: chi(p, x), bounds=(-13, -8), method="bounded")
    res = go - np.log10(gb*(1 + (10**b.x/gb)**(2*p)))
    print(f"   p = {p:4}: {lab:38s} best a_c = {10**b.x:.2e}  chi2 = {b.fun:7.0f}  scatter {np.std(res):.3f} dex  trend {np.polyfit(np.log10(gb), res, 1)[0]:+.2f}")
b = minimize(lambda q: chi(q[0], q[1]), [0.3, -10], method="Nelder-Mead")
print(f"   free fit: p = {b.x[0]:.3f}, a_c = {10**b.x[1]:.2e}, chi2 = {b.fun:.0f}")
ref = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
pred = np.log10(gb*ref(gb/1.15e-10)); print(f"   (reference: the slack rule used before, chi2 = {np.sum(((go-pred)/np.hypot(el,0.1))**2):.0f})")
