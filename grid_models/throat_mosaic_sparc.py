# SPARC + Cassini for the mosaic-derived switch-over mu(x) (one free scale), vs McGaugh; low-pull slope of the RAR
import numpy as np, json, glob
from scipy.optimize import minimize_scalar
from scipy.interpolate import PchipInterpolator
src = open("taut_refill.py").read()
exec(src.split("mus = {")[0]); exec("def make_nu" + src.split("def make_nu")[1].split("nus = {k")[0])
Y = np.logspace(-8, 8, 4001); exec("def Q2" + src.split("def Q2")[1].split('print(f"\\nCassini')[0])
GM = 1.32712e20; gext = 2.32e-10; kpc = 3.0857e19
M = json.load(open("throat_from_mosaic.json"))
nus = {"McGaugh 2016": lambda y: 1/(1-np.exp(-np.sqrt(y)))}
for rule, d in M.items():
    t, K = np.array(d["t"]), np.array(d["K"]); x = t - d["tc"]; m = x >= 0
    xs, Ks = x[m], np.maximum.accumulate(K[m]); xs, idx = np.unique(xs, return_index=True); Ks = Ks[idx]
    f = PchipInterpolator(xs, Ks)
    mu = (lambda f, xmax: (lambda z: np.where(z < xmax, np.clip(f(np.minimum(z, xmax)), 1e-12, 1), 1.0)))(f, xs[-1])
    nus[f"mosaic, {rule}"] = make_nu(mu)
GO, GB, EL = [], [], []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    dd = np.loadtxt(fn, comments="#", ndmin=2)
    if np.any(dd[:, 5] != 0): continue
    r, V, eV, Vg, Vd, Vb, _, _ = dd.T; gb = (Vg*abs(Vg)+0.6*Vd**2)*1e6/(r*kpc); go = V**2*1e6/(r*kpc)
    ok = (V > 0)&(eV/V < 0.1)&(gb > 0)&(r > 0); GO += list(go[ok]); GB += list(gb[ok]); EL += list(2*eV[ok]/V[ok]/np.log(10))
go, gb, el = map(np.array, (GO, GB, EL)); lo = gb < 3e-11
print(f"bulge-free SPARC, M/L 0.6: {len(go)} pts. Measured low-pull slope d log g_obs / d log g_bar (g_bar<3e-11): {np.polyfit(np.log10(gb[lo]), np.log10(go[lo]), 1)[0]:.2f}  (square-root law = 0.50)")
ref = None
for lab, nu in nus.items():
    chi = lambda a: np.sum(((np.log10(go)-np.log10(gb*nu(gb/a)))/np.hypot(el, 0.1))**2)
    b = minimize_scalar(lambda z: chi(10**z), bounds=(-12, -8), method="bounded"); a = 10**b.x
    ref = b.fun if ref is None else ref
    pred = np.log10(gb*nu(gb/a)); sl = np.polyfit(np.log10(gb[lo]), pred[lo], 1)[0]
    print(f"   {lab:18s} scale {a:.2e}  dchi2 vs McGaugh {b.fun-ref:+8.1f}  predicted low-pull slope {sl:.2f}  Cassini Q2 {Q2(nu, a)*1e27:6.1f}e-27")
