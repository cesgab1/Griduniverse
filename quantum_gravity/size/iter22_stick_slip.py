"""
ITERATION 22 (Coalesce: stretchable things have built-in tension, like rubber). Stick-slip grid:
 each cell's strain grows at the stretch rate (d eps/dt = H); when it reaches the yield threshold eps_c a new cell is added
 and the strain drops by eps_c. Optional relief delay tau (the new cell appears tau after the threshold is reached).
 Energy density of the leftover strain: rho = K <eps^2>/2, K = grid stiffness (~ Planck density).
A. Simulation (scaled: eps_c = 1e-3, 20000 cells, random phases) through a matter -> Lambda-like history of H:
   is <eps^2> constant (frozen size) or does it follow H?
B. With a relief delay tau: <eps^2> = eps_c^2/3 + eps_c H tau + (H tau)^2 (overshoot): which way does the size move?
C. What threshold gives the observed size, and what would set it?
D. Data: what history does this give, compared with Claim 1 and the jostled toy?
"""
import numpy as np
rng = np.random.default_rng(22); out = ["ITERATION 22: stick-slip grid (built-in tension with a yield threshold)", ""]
Om = 0.31; lna = np.linspace(np.log(1/30), 0, 400); a = np.exp(lna); H = np.sqrt(Om*a**-3 + 1 - Om)
t = np.r_[0, np.cumsum(np.diff(lna)/H[1:])]
def run(eps_c, tau, ncell=20000):
    eps = rng.uniform(0, eps_c, ncell); pend = np.full(ncell, np.inf); res = []
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]; nsub = max(1, int(np.ceil(H[i]*dt/(eps_c/20)))); h = dt/nsub
        for _ in range(nsub):
            eps += H[i]*h
            hit = (eps >= eps_c) & np.isinf(pend); pend[hit] = tau
            pend[~np.isinf(pend)] -= h
            done = pend <= 0; eps[done] -= eps_c; pend[done] = np.inf
        res.append(np.mean(eps**2))
    return np.array(res)
eps_c = 1e-3
for tau in (0.0, 2e-6, 5e-6, 1e-3):
    r = run(eps_c, tau); rel = r/(eps_c**2/3)
    out.append(f"A/B. relief delay tau = {tau:6.0e} (Hubble times): <eps^2>/(eps_c^2/3) at z = 29 / 3 / 1 / 0 = "
               + " / ".join(f"{rel[np.argmin(abs(1/a[1:] - 1 - zz))]:.3f}" for zz in (29, 3, 1, 0)))
out += ["   -> prompt relief: the leftover strain is CONSTANT (frozen size, exactly like a cosmological constant, w = -1).",
        "   -> short delay (tau << eps_c/H, here eps_c/H ~ 1e-5 early): overshoot ∝ H tau that SHRINKS as the expansion slows",
        "      (larger early, smaller late): w > -1 at all times, never crossing.",
        "   -> long delay (tau > eps_c/H): relief cannot keep up, strain piles up without limit (runaway) - not a viable grid.", ""]
rhoP, rho_obs = 5.16e96, 5.91e-27; K = rhoP
eps_need = np.sqrt(3*2*rho_obs/K); lP = 1.616e-35; LH = 2.998e8/(67.4e3/3.0857e22)
out.append(f"C. size needs eps_c = {eps_need:.1e}. If one new cell relieves strain along a whole row of N cells, eps_c = 1/N:")
out.append(f"   N = {1/eps_need:.1e} cells = a row {1/eps_need*lP:.1e} m long; today's Hubble length = {LH:.1e} m (ratio {1/eps_need*lP/LH:.1f})")
out.append("   -> the threshold is 'one cell of stretch per horizon-long row'. If that row follows the CURRENT horizon, the size tracks")
out.append("      H^2 (iteration 21's failure); if it was fixed once, its length must equal today's horizon by coincidence ('why now').")
out += ["", "D. Data: prompt stick-slip = a pure cosmological constant -> Delta chi2 = 0 vs Lambda by construction, i.e. it LOSES Claim 1's",
        "   advantage (-5.4 / -7.2 / -6.8). As a positive baseline under the jostled toy it is disfavoured (iteration 13: data want the",
        "   jostled fraction f >= 1, baseline <= 0). With a relief delay it is w > -1 throughout (thawing-like), against the crossing",
        "   the data prefer.",
        "", "Lesson: built-in tension with a yield threshold gives a FROZEN size naturally (Coalesce's point is right), but its natural",
        "behaviour is a pure constant; it does not carry Claim 1's shape, and the threshold's smallness is the old problem in a new",
        "form ('one cell of stretch per horizon-long row')."]
txt = "\n".join(out); print(txt); open("iter22_stick_slip.txt", "w").write(txt + "\n")
