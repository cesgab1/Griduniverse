"""
Slack from drained, thinned, lengthened links -> does it give the galaxy law?

Picture: where the fluid has drained, layers collapse, links thin and (Poisson) lengthen, so each link carries a
little excess length s ("slack"): it hangs loose until the local stretch exceeds s, then pulls linearly.
A patch of grid = many links with a random spread of slack p(s) (random mosaic -> random s).
Restoring stress of the patch at local stretch e:   sigma(e) = k * Int_0^e (e - s) p(s) ds      (fibre recruitment)
Gravity: stretch e <-> local pull g (slope of the dimple); the mass sets sigma = g_N (flux of the source).
Nothing is fitted except the slack scale; the SHAPE follows from p(s).

Deep regime (e << typical s):  sigma ~ p(0) e^2 / 2   ->  g = sqrt(2 g_N / p(0)) = sqrt(a0 g_N),  a0 = 2/p(0)
Strong regime (e >> all s):    sigma ~ e - <s>        ->  g = g_N + <s>   (a leftover constant pull)
"""
import numpy as np, glob, os
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erf

# slack distributions, each written with a0 = 2/p(0) as the single scale
def sig_uniform(e, a0):          # p = 1/sm on [0, sm], sm = a0/2
    sm = a0/2; e = np.asarray(e, float)
    return np.where(e < sm, e**2/(2*sm), e - sm/2)
def sig_expo(e, a0):             # p = exp(-s/s0)/s0, s0 = a0/2
    s0 = a0/2; return e - s0*(1 - np.exp(-e/s0))
def sig_halfgauss(e, a0):        # p = sqrt(2/pi)/w exp(-s^2/2w^2), a0 = sqrt(2 pi) w
    w = a0/np.sqrt(2*np.pi); x = e/(np.sqrt(2)*w)
    # Int_0^e (e-s) p ds
    return e*erf(x) - w*np.sqrt(2/np.pi)*(1 - np.exp(-x**2))
models = {"uniform slack": (sig_uniform, 0.25), "exponential slack": (sig_expo, 0.5),
          "half-Gaussian slack": (sig_halfgauss, 1/np.pi)}   # offset/a0 = <s> p(0)/2

def g_of_gN(gN, a0, sig):
    out = np.empty_like(gN)
    for i, y in enumerate(gN):
        out[i] = brentq(lambda g: sig(g, a0) - y, 1e-20, y + 10*a0 + 1e-12)
    return out

# ---------------- SPARC radial acceleration relation ----------------
Yd, Yb, kpc = 0.5, 0.7, 3.0857e19
GO, GB, EL = [], [], []
for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(f, comments="#");  d = d[None] if d.ndim == 1 else d
    r, V, eV, Vg, Vd, Vb, SBd, SBb = d.T
    gb = (Vg*np.abs(Vg) + Yd*Vd**2 + Yb*Vb**2)*1e6/(r*kpc); go = V**2*1e6/(r*kpc)
    ok = (V > 0) & (eV/V < 0.1) & (gb > 0) & (r > 0)
    GO += list(go[ok]); GB += list(gb[ok]); EL += list(2*eV[ok]/V[ok]/np.log(10))
go, gb, el = map(np.array, (GO, GB, EL)); lgo, lgb = np.log10(go), np.log10(gb)
print(f"SPARC: {len(go)} points\n")

def score(pred):
    res = lgo - np.log10(pred); return np.sum((res/np.hypot(el, 0.1))**2), np.std(res), res
refs = {"McGaugh 2016 empirical": lambda gN, a0: gN/(1 - np.exp(-np.sqrt(gN/a0)))}
print(f"{'law':28s} {'best a0':>10s} {'chi2':>8s} {'scatter':>8s} {'trend':>7s}")
results = {}
for lab, fn in list(refs.items()) + [(k, (lambda s: (lambda gN, a0: g_of_gN(gN, a0, s)))(v[0])) for k, v in models.items()]:
    # coarse-to-fine a0 search
    grid = np.linspace(-10.4, -9.4, 41)
    chis = [score(fn(gb, 10**x))[0] for x in grid]
    x0 = grid[int(np.argmin(chis))]
    b = minimize_scalar(lambda x: score(fn(gb, 10**x))[0], bounds=(x0-0.03, x0+0.03), method="bounded")
    a0 = 10**b.x; chi, sc, res = score(fn(gb, a0)); tr = np.polyfit(lgb, res, 1)[0]
    results[lab] = a0
    print(f"{lab:28s} {a0:10.3e} {chi:8.1f} {sc:8.3f} {tr:+7.3f}")

# residuals in bins of pull, to see WHERE each law misses
bins = np.arange(-12, -8.4, 0.5)
print("\nmedian residual (dex) by pull bin   log g_bar:", " ".join(f"{b+0.25:6.2f}" for b in bins[:-1]))
for lab in results:
    fn = refs.get(lab) or (lambda gN, a0, s=models[lab][0]: g_of_gN(gN, a0, s))
    res = lgo - np.log10(fn(gb, results[lab]))
    row = [np.median(res[(lgb >= lo) & (lgb < lo+0.5)]) if np.any((lgb >= lo) & (lgb < lo+0.5)) else np.nan for lo in bins[:-1]]
    print(f"   {lab:28s}", " ".join(f"{v:+6.3f}" for v in row))

# ---------------- a0 from the grid (carried from the Verlinde elastic derivation) ----------------
c, H0 = 2.998e8, 67.5e3/3.0857e22
print(f"\ncH0/6 = {c*H0/6:.3e} m/s^2.  In stretch units e = g/(cH0): a0/(cH0) = 1/6 -> p(0) = 12,")
print(f"   i.e. the loosest links carry ~{1/12*100:.0f}% excess length (uniform) -- needs to come out of the drainage/Poisson geometry; NOT derived yet.")

# ---------------- Solar System: the leftover constant pull ----------------
GM, a_sat, e_sat = 1.32712e20, 1.4335e12, 0.0565
P_yr = 2*np.pi*np.sqrt(a_sat**3/GM)/3.156e7
print("\nSolar System check (Saturn, g_N = %.2e m/s^2):" % (GM/a_sat**2))
for lab, (s, frac) in models.items():
    a0 = results[lab]; dg = frac*a0
    # constant radial perturbing acceleration A -> perihelion shift per orbit = 2 pi A sqrt(1-e^2) a^2 / GM (retrograde for inward A)
    dw = 2*np.pi*dg*np.sqrt(1-e_sat**2)*a_sat**2/GM
    print(f"   {lab:22s} leftover pull {dg:.1e} m/s^2 -> Saturn perihelion drift {dw/P_yr*206264.8e3:7.1f} mas/yr")
print("   ephemeris bounds on Saturn's anomalous perihelion drift are at the ~1 mas/yr level (INPOP/EPM).")
# theorem: for any slack spread that is densest at zero slack, <s> >= 1/(2 p(0)) -> leftover >= a0/4
print("   Any slack spread densest at zero slack gives leftover >= a0/4 (uniform is the minimum).")
