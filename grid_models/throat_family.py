"""
With the drained grid at its rigidity point at zero pull (sealed fraction f0 = f_r) and stiffness linear in (f - f_r)
(measured, drained_stiffness.txt), the refill rule gives mu(x) = C(x): the cumulative spread of throat strengths
(how hard it is for fluid to leave a gap). Deep MOND for any C with finite density at zero.
Question: is there ANY spread with a hard upper end that satisfies both SPARC and Cassini?
Family: C(x) = 1-(1-x/n)^n for x<n, 1 above  (a0 = 1/C'(0); all fluid resealed at pull n*a0).
n=1 uniform ... n->inf exponential.  Not a fit: a scan to map the allowed region.
"""
import numpy as np, glob
from scipy.optimize import minimize_scalar
src = open("taut_refill.py").read()
exec(src.split("mus = {")[0])
exec("def make_nu" + src.split("def make_nu")[1].split("nus = {k")[0])
Y = np.logspace(-8, 8, 4001)
exec("def Q2" + src.split("def Q2")[1].split('print(f"\\nCassini')[0])
GM = 1.32712e20; gext = 2.32e-10; kpc = 3.0857e19; c, H0 = 2.998e8, 67.5e3/3.0857e22; a0g = c*H0/6
fam = {f"n={n}": (lambda n: (lambda x: np.where(x < n, 1-(1-np.minimum(x, n)/n)**n, 1.0)))(n) for n in (1, 1.5, 2, 3, 5)}
fam["n=inf (exponential)"] = lambda x: 1-np.exp(-x)
nus = {k: make_nu(m) for k, m in fam.items()}
nus["McGaugh 2016 (reference)"] = lambda y: 1/(1-np.exp(-np.sqrt(y)))
def data(Yd, nobulge):
    GO, GB, EL = [], [], []
    for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
        d = np.loadtxt(f, comments="#", ndmin=2)
        if nobulge and np.any(d[:, 5] != 0): continue
        r, V, eV, Vg, Vd, Vb, _, _ = d.T
        gb = (Vg*abs(Vg)+Yd*Vd**2+1.4*Yd*Vb**2)*1e6/(r*kpc); go = V**2*1e6/(r*kpc)
        ok = (V > 0)&(eV/V < 0.1)&(gb > 0)&(r > 0); GO += list(go[ok]); GB += list(gb[ok]); EL += list(2*eV[ok]/V[ok]/np.log(10))
    return map(np.array, (GO, GB, EL))
sets = {"all, M/L .5": data(.5, 0), "all, M/L .6": data(.6, 0), "no-bulge, M/L .5": data(.5, 1), "no-bulge, M/L .6": data(.6, 1)}
sets = {k: tuple(v) for k, v in sets.items()}
print("chi2 is relative to the McGaugh curve in the same sample (negative = better); Q2 at that sample's best a0 (bound 3 +- 3)")
hdr = "".join(f"{k:>24s}" for k in sets); print(f"{'':26s}{hdr}\n{'':26s}" + "   dchi2  a0e10   Q2e27 "*len(sets))
ref = {}
for lab in ["McGaugh 2016 (reference)"] + list(fam):
    nu = nus[lab]; line = f"{lab:26s}"
    for sk, (go, gb, el) in sets.items():
        chi = lambda a0: np.sum(((np.log10(go)-np.log10(gb*nu(gb/a0)))/np.hypot(el, 0.1))**2)
        b = minimize_scalar(lambda x: chi(10**x), bounds=(-10.4, -9.4), method="bounded")
        if lab.startswith("McGaugh"): ref[sk] = b.fun
        line += f"  {b.fun-ref[sk]:+7.1f} {10**b.x*1e10:5.2f} {Q2(nu, 10**b.x)*1e27:7.1f} "
    print(line, flush=True)
