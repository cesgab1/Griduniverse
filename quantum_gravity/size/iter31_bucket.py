"""
ITERATION 31: Coalesce's 'time adds it' as a BUCKET fed by the SAME random kicks as the jostled tension (the 'sponge').
Toy (iteration 1/4): tension tau gets kicks with variance rate 1/a per unit time; the sponge forgets at rate kappa H (kappa = 3,
Hubble friction); energy = A * variance (quadratic).  Bucket: the same kicks, NO forgetting:
      dX_J/dlna = -2 kappa X_J + 1/(aH)        (sponge, as before)
      dX_B/dlna =                1/(aH)        (bucket: X_B = conformal time since the start)
      rho_DE = A (X_J + X_B),  A fixed by today's dark energy.   NO new free parameter (the split is fixed by the physics).
Expectations written BEFORE fitting:
 E1 sign: energy is quadratic, so the bucket is POSITIVE. The slightly negative accumulated part the data mildly like
    (iteration 13) is unreachable in this toy -> FAIL on sign, by construction.
 E2 size: today X_B/X_J ~ (conformal age ~3.2/H0) / (1/(6 H0)) ~ 19: the bucket dominates (~95% of dark energy).
 E3 history: dominated by X_B ~ conformal time, growing -> w ~ -1.1 today (phantom now), the opposite of the data's w > -1.
    Predicted: Delta chi2 vs Lambda POSITIVE (worse) for all three supernova sets.
"""
import numpy as np, os, sys
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
here = os.path.dirname(os.path.abspath(__file__))
tg = os.path.join(here, "..", "toy_grid", "iter1_memory.py")
src = open(tg).read().split("# ---------- Part A")[0]
g = {"__name__": "x", "__file__": os.path.abspath(tg)}; cwd = os.getcwd(); os.chdir(os.path.dirname(tg)); exec(src, g); os.chdir(cwd)
fl = g["g"]; ZG, W_R = fl["ZG"], fl["W_R"]
MODE = os.environ.get("BUCKET", "on")
kappa = 3.0
def shape(Om, h, a_grid):
    Or = W_R/h**2; OL = 1 - Om - Or; ai = 1/31.
    etai = quad(lambda a: 1/(a*a*np.sqrt(Om*a**-3 + Or*a**-4)), 1e-9, ai)[0]   # conformal time before z = 30 (DE negligible)
    def run(A):
        def rhs(lna, y):
            a = np.exp(lna); XJ, XB = max(y[0], 1e-300), y[1]
            rde = A*(XJ + (XB if MODE == "on" else 0.0))
            H = np.sqrt(Om*a**-3 + Or*a**-4 + rde)
            return [-2*kappa*XJ + 1/(a*H), 1/(a*H)]
        Hi = np.sqrt(Om*ai**-3 + Or*ai**-4)
        return solve_ivp(rhs, [np.log(ai), 0.0], [1/(2*kappa*ai*Hi), etai], t_eval=np.log(a_grid[::-1]), rtol=1e-8, atol=1e-14)
    def today(A):
        s = run(A); return A*(s.y[0][-1] + (s.y[1][-1] if MODE == "on" else 0)) - OL
    A = brentq(today, 1e-8, 1e6, xtol=1e-12)
    s = run(A); XJ, XB = s.y[0][::-1], s.y[1][::-1]
    return A*(XJ + (XB if MODE == "on" else 0)), XB[0]/XJ[0]
def E_b(model, Om, h, extra):
    Or = W_R/h**2; zp = 1 + ZG; zi = ZG[ZG <= 30]; a = 1/(1 + zi)
    r, _ = shape(Om, h, a); rde = np.zeros_like(ZG); rde[:len(zi)] = r
    # beyond z = 30: bucket keeps its (tiny) value; negligible
    return np.sqrt(Om*zp**3 + Or*zp**4 + rde)
if os.environ.get("REPORT"):
    a = 1/(1 + np.linspace(0, 3, 301)); r, ratio = shape(0.31, 0.68, a)
    z = 1/a - 1; w = -1 - np.gradient(np.log(r), np.log(a))/3
    print(f"bucket/sponge today = {ratio:.1f}; w0 = {w[0]:+.3f}; w(z=0.5) = {np.interp(0.5, z, w):+.3f}; w(z=1) = {np.interp(1, z, w):+.3f}; w(z=2) = {np.interp(2, z, w):+.3f}")
    sys.exit()
L = fl["best"]("LCDM", [[0.31, 0.68, 0.0224]])
fl["E_of_z"] = E_b
r = fl["best"]("B", [[0.31, 0.68, 0.0224]])
print(f"RESULT {os.environ.get('SNSET', 'PANTHEON')} bucket={MODE}: dchi2 vs LCDM {r.fun - L.fun:+.2f}")
