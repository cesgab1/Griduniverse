"""
ITERATION 1 of the grid toy: does the grid tension have to follow the stretch rate INSTANTLY, or may it lag?

Toy rules (stated, not derived):
 - Space is filled by fixed Planck cells (added by division as space grows; daughters inherit the parent's tension).
 - Comoving domains (size L comoving, a*L physical) carry an intensive tension tau.
 - Kicks: flows crossing a domain hit its walls at rate c/(a L): each hit adds ±eps to tau (random sign).
 - Memory: tension relaxes at rate kappa*H (memory time 1/(kappa H), a fraction of the Hubble time).
 => d<tau^2>/dt = -2 kappa H <tau^2> + eps^2 c/(a L).
 - Energy density of the tension: LINEAR rho ∝ <|tau|> ∝ sqrt(<tau^2>)  (string-like)   -> instant limit: Claim 1 (adot^-1/2)
                                   QUADRATIC rho ∝ <tau^2>                (spring-like)  -> instant limit: adot^-1   (beta = 1)
Cell division with inheritance copies tau, so it does not change the distribution of tau: the grid reduces exactly to one
stochastic variable per domain. Part A checks that with a direct Monte Carlo (many domains). Part B solves the
self-consistent expansion (H depends on rho_DE) for each memory kappa and fits DESI DR2 + Planck priors + supernovae.
"""
import numpy as np, os, sys, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
FIT = "/home/claude/griduniverse/cosmology_fits/fit_law.py"
os.environ["ONLY"] = "NONE"
g = {"__name__": "fitlaw", "__file__": FIT}
try: exec(open(FIT).read(), g)
except SystemExit: pass
ZG, W_R = g["ZG"], g["W_R"]

def rho_hist(Om, h, kappa, energy, a_grid):
    """self-consistent rho_DE(a) in units of today's critical density; rho_DE(1) = 1 - Om - Or."""
    Or = W_R/h**2; OL = 1 - Om - Or; ai = 1/31.
    def run(A):
        def rhs(lna, y):
            a = np.exp(lna); X = max(y[0], 1e-300)
            rde = A*np.sqrt(X) if energy == "lin" else A*X
            H = np.sqrt(Om*a**-3 + Or*a**-4 + rde)
            if kappa == np.inf: return [0.0]
            return [-2*kappa*X + 1/(a*H)]
        if kappa == np.inf: return None
        Hi = np.sqrt(Om*ai**-3 + Or*ai**-4); X0 = 1/(2*kappa*ai*Hi)     # start in balance (dark energy negligible at z = 30)
        s = solve_ivp(rhs, [np.log(ai), 0.0], [X0], t_eval=np.log(a_grid[::-1]) if a_grid is not None else [0.0], rtol=1e-8, atol=1e-14)
        return s
    def today(A):
        s = run(A); X = s.y[0][-1]; return (A*np.sqrt(X) if energy == "lin" else A*X) - OL
    # bracket A
    lo, hi = 1e-6, 1e6
    A = brentq(today, lo, hi, xtol=1e-10, rtol=1e-10, maxiter=200)
    s = run(A); X = s.y[0][::-1]
    return A*np.sqrt(X) if energy == "lin" else A*X

def make_E(kappa, energy):
    def E_of_z(model, Om, h, extra):
        Or = W_R/h**2; zp = 1 + ZG
        zi = ZG[ZG <= 30]; a = 1/(1 + zi)
        rde = np.zeros_like(ZG); rde[:len(zi)] = rho_hist(Om, h, kappa, energy, a)
        return np.sqrt(Om*zp**3 + Or*zp**4 + rde)
    return E_of_z

# ---------- Part A: Monte Carlo check of the reduction (5000 domains, real kicks, matter + Lambda background) ----------
def part_A(kappa=1.0, n=5000, steps=6000):
    rng = np.random.default_rng(1); Om = 0.31; lna = np.linspace(np.log(1/31), 0, steps); d = lna[1] - lna[0]
    H = lambda a: np.sqrt(Om*a**-3 + 1 - Om); a0 = np.exp(lna[0])
    tau = rng.normal(0, np.sqrt(1/(2*kappa*a0*H(a0))), n); X = 1/(2*kappa*a0*H(a0)); res = []
    for k in range(steps - 1):
        a = np.exp(lna[k]); dt = d/H(a)                                   # time step (units 1/H0)
        tau += -kappa*H(a)*tau*dt + np.sqrt(dt/a)*rng.standard_normal(n)  # Euler-Maruyama; kick variance rate 1/a
        X += (-2*kappa*X + 1/(a*H(a)))*d
        if k % 1500 == 0 or k == steps - 2: res.append((1/a - 1, np.mean(tau**2), X))
    return res

out = ["ITERATION 1: grid tension with finite memory (kappa = memory rate / H)", "", "Part A: Monte Carlo of 5000 domains vs the one-variable equation (kappa = 1):"]
for z, mc, ode in part_A(): out.append(f"  z = {z:6.2f}: Monte Carlo <tau^2> = {mc:9.4f}   equation {ode:9.4f}   ratio {mc/ode:.3f}")
print("\n".join(out)); sys.stdout.flush()

# ---------- Part B: fits ----------
SN = os.environ.get("SNSET", "PANTHEON")
res = {}
L = g["best"]("LCDM", [[0.31, 0.68, 0.0224], [0.30, 0.69, 0.0223]]); res["LCDM"] = float(L.fun)
for energy in ("lin", "quad"):
    for kappa in (0.25, 0.5, 1.0, 2.0, 4.0, 16.0):
        g["E_of_z"] = make_E(kappa, energy)
        r = g["best"]("MEM", [[0.31, 0.68, 0.0224]])
        res[f"{energy} kappa={kappa}"] = (float(r.fun), [float(x) for x in r.x])
        print(f"{SN} {energy:4s} kappa {kappa:5.2f}: dchi2 vs LCDM {r.fun - L.fun:+6.2f}  params {np.round(r.x, 4)}"); sys.stdout.flush()
json.dump(res, open(os.path.join(HERE, f"iter1_fits_{SN}.json"), "w"), indent=1)
open(os.path.join(HERE, "iter1_partA.txt"), "w").write("\n".join(out) + "\n")
