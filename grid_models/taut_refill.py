"""
Refill rule (no new inputs): a link that is stretched past its slack goes taut, seals the gap between layers,
and the fluid flows back in -> its Poisson lengthening (slack) is undone. Loose links stay drained.
Patch stress:  sigma(e) = k * e * F(e),  F = cumulative slack distribution (fraction of links already taut).
=> MOND interpolating function mu(x) = F(x) exactly.  Deep: sigma ~ p(0) e^2 -> sqrt law, a0 = 1/p(0).
Strong pull: all links taut -> sigma = k e -> pure Newton, NO leftover pull.
Exponential slack (random-mosaic chord lengths): mu = 1 - exp(-x).
Tests: SPARC RAR; Cassini quadrupole Q2 (QUMOND) with Galactic field g_ext = 2.32e-10 (Gaia EDR3).
"""
import numpy as np, glob
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erf
from scipy.interpolate import interp1d

mus = {"exponential slack  mu=1-exp(-x)": lambda x: 1-np.exp(-x),
       "half-Gaussian slack mu=erf(x/sqrt2)": lambda x: erf(x/np.sqrt(2)),
       "uniform slack  mu=min(x,1)": lambda x: np.minimum(x, 1.0)}
# nu(y): g = nu * gN, with gN = g mu(g)   (units of the scale)
Y = np.logspace(-8, 8, 4001)
def make_nu(mu):
    x = np.array([brentq(lambda x: x*mu(x) - y, 1e-12, y + 50) for y in Y])
    lnnu = interp1d(np.log(Y), np.log(x/Y), fill_value="extrapolate")
    return lambda y: np.exp(lnnu(np.log(np.clip(y, 1e-8, 1e8))))
nus = {k: make_nu(m) for k, m in mus.items()}
nus["McGaugh 2016 empirical (reference)"] = lambda y: 1/(1-np.exp(-np.sqrt(y)))

# ---------------- SPARC ----------------
kpc = 3.0857e19
D = [np.loadtxt(f, comments="#", ndmin=2) for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat"))]
def load(Yd):
    GO, GB, EL = [], [], []
    for d in D:
        r, V, eV, Vg, Vd, Vb, _, _ = d.T
        gb = (Vg*abs(Vg) + Yd*Vd**2 + 1.4*Yd*Vb**2)*1e6/(r*kpc); go = V**2*1e6/(r*kpc)
        ok = (V > 0) & (eV/V < 0.1) & (gb > 0) & (r > 0)
        GO += list(go[ok]); GB += list(gb[ok]); EL += list(2*eV[ok]/V[ok]/np.log(10))
    return map(np.array, (GO, GB, EL))
c, H0 = 2.998e8, 67.5e3/3.0857e22; a0g = c*H0/6
best = {}
for Yd in (0.5, 0.6):
    go, gb, el = load(Yd)
    print(f"\nSPARC, stellar mass-to-light {Yd} ({len(go)} pts)   law | best a0 | chi2 | scatter | chi2 at a0=cH0/6")
    for lab, nu in nus.items():
        chi = lambda a0: np.sum(((np.log10(go) - np.log10(gb*nu(gb/a0)))/np.hypot(el, 0.1))**2)
        b = minimize_scalar(lambda x: chi(10**x), bounds=(-10.4, -9.4), method="bounded")
        a0 = 10**b.x; res = np.log10(go) - np.log10(gb*nu(gb/a0))
        best[(lab, Yd)] = a0
        print(f"   {lab:38s} {a0:.3e}  {b.fun:7.1f}  {np.std(res):.3f}  {chi(a0g):7.1f}")

# ---------------- Solar System ----------------
GM = 1.32712e20; gext = 2.32e-10
def Q2(nu, a0):
    rM = np.sqrt(GM/a0)
    ge = gext/a0
    # Newtonian external field (dimensionless): solve ge = nu(yN) yN
    eta = brentq(lambda y: nu(y)*y - ge, 1e-6, 1e3)
    lr = np.linspace(np.log(1e-3), np.log(1e3), 3000); r = np.exp(lr)
    ct = np.linspace(-1, 1, 801)
    R, C = np.meshgrid(r, ct, indexing="ij")
    # gN = -rhat/r^2 + eta zhat ; components in (r, theta): 
    gr = -1/R**2 + eta*C; gt = -eta*np.sqrt(1-C**2)
    y = np.hypot(gr, gt)
    f = nu(y) - 1
    # S = div[(nu-1) gN] = gN . grad(nu) ; compute via finite differences of f
    dfdr = np.gradient(f, r, axis=0)
    th = np.arccos(C); dfdth = np.gradient(f, ct, axis=1)*(-np.sqrt(1-C**2))  # d/dtheta = -sin * d/dcos
    S = gr*dfdr + gt*dfdth/R
    # Phi_Q = -Q2/2 x_i x_j (e_i e_j - d_ij/3) -> d2Phi/dz2 = -2Q2/3 ;  lap(dPhi) = -S (g=-grad Phi)
    # dPhi = (1/4pi) Int S/|x-x'| -> d2/dz2 at 0 = (1/4pi) Int S * 2 P2(cos)/r^3 d3x
    P2 = 0.5*(3*C**2 - 1)
    integrand = S*2*P2/R**3 * R**2 * 2*np.pi / (4*np.pi)
    d2 = np.trapezoid(np.trapezoid(integrand, ct, axis=1), r)
    return -1.5*d2 * a0/rM
print(f"\nCassini quadrupole: bound Q2 = (3 +- 3)e-27 s^-2 ; Galactic field {gext:.2e} m/s^2")
for lab, nu in nus.items():
    for Yd in (0.5, 0.6):
        a0 = best[(lab, Yd)]
        print(f"   {lab:38s} a0={a0:.2e} (M/L {Yd}): Q2 = {Q2(nu, a0)*1e27:7.2f}e-27")
    print(f"   {'':38s} a0=cH0/6            : Q2 = {Q2(nu, a0g)*1e27:7.2f}e-27")
