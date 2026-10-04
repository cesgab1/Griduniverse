"""
ITERATION 50 (pre-registered in PREREG_50.md): VCDM-type linear growth for Claim 1 vs the same background with the ordinary
constraint, and vs LCDM. Units H0 = 1, rho in units of 3 M_P^2 (H^2 = rho_total). Matter-dominated start (z = 200),
radiation kept in H but not perturbed (same in all runs; ratios are what matter).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
# REVISION after the first run (recorded, output kept in iter50_first_run_with_radiation.txt): P1 failed (|F-1| = 0.075 on
# LCDM) because radiation entered dH/dt but was not a perturbed fluid in F's denominator. Fix: drop radiation everywhere in this
# calculation (start z = 200, where it is ~6% of matter; the same in all runs).
Om, Or = 0.315, 0.0; Ode = 1 - Om - Or; hubble_mpc = 2997.9   # c/H0 in Mpc/h
def E_claim1(a):
    g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - Ode*(a*np.exp(le))**-0.5
    return np.exp(brentq(g, -20, 60, xtol=1e-14, maxiter=500))
def E_lcdm(a): return np.sqrt(Om*a**-3 + Or*a**-4 + Ode)
lna_grid = np.linspace(np.log(1/201), 0, 3001)
def tables(Efun):
    E = np.array([Efun(np.exp(x)) for x in lna_grid]); dlnE = np.gradient(np.log(E), lna_grid)
    return E, dlnE
TAB = {"claim1": tables(E_claim1), "lcdm": tables(E_lcdm)}
def run(k_hmpc, bg, vcdm):
    k = k_hmpc*hubble_mpc; E_t, dlnE_t = TAB[bg]
    def rhs(x, y):
        Phi, d, th = y; a = np.exp(x); E = np.interp(x, lna_grid, E_t); dHdt = E*E*np.interp(x, lna_grid, dlnE_t)
        cH = a*E; rm = Om*a**-3
        F = (k**2 - 3*a*a*dHdt)/(k**2 + 4.5*a*a*rm) if vcdm else 1.0
        dPhi = -Phi + 1.5*a*a*rm*th*F/(k**2*cH)
        return [dPhi, -th/cH + 3*dPhi, -th + k**2*Phi/cH]
    x0 = lna_grid[0]; a0 = np.exp(x0); cH0 = a0*TAB[bg][0][0]; Phi0 = 1.0
    d0 = -2*Phi0 - (2/3)*k**2*Phi0/cH0**2; th0 = (2/3)*k**2*Phi0/cH0
    s = solve_ivp(rhs, [x0, 0], [Phi0, d0, th0], rtol=1e-10, atol=1e-14, dense_output=True)
    Phi, d, th = s.y[:, -1]; cHt = TAB[bg][0][-1]
    return d + 3*cHt*th/k**2, Phi, s                   # comoving density contrast, potential today
out = ["ITERATION 50: horizon-scale growth with the frame field (expectations pre-registered)", ""]
# P1: F on LCDM background
E_t, dl_t = TAB["lcdm"]; a = np.exp(lna_grid); k = 0.001*hubble_mpc
F = (k**2 - 3*a*a*E_t**2*dl_t)/(k**2 + 4.5*a*a*Om*a**-3)
out.append(f"P1 F on the LCDM background: max |F - 1| = {np.max(np.abs(F - 1)):.1e} (should be ~0; radiation term and finite differences)")
out += ["", "   k [h/Mpc]   A/B density today   A/B potential today   (A = our model, B = same background, ordinary constraint)"]
rows = []
for kk in (1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1):
    DA, PA, _ = run(kk, "claim1", True); DB, PB, _ = run(kk, "claim1", False); DC, PC, _ = run(kk, "lcdm", False)
    rows.append((kk, DA/DB, PA/PB, DB/DC, PB/PC))
    out.append(f"   {kk:9.0e}   {DA/DB:.6f}            {PA/PB:.6f}")
out += ["", "   for comparison, background effect alone (B vs LCDM): density ratio at k = 0.01: "
        f"{[r for r in rows if r[0] == 1e-2][0][3]:.4f}, potential ratio {[r for r in rows if r[0] == 1e-2][0][4]:.4f}"]
# ISW proxy: potential decay since z = 10, largest scales
for kk in (2e-4, 1e-3, 1e-2):
    res = {}
    for lab, bg, v in (("A", "claim1", True), ("B", "claim1", False), ("C", "lcdm", False)):
        _, P, s = run(kk, bg, v); x10 = np.log(1/11); res[lab] = 1 - P/s.sol(x10)[0]
    out.append(f"   potential decay z=10 -> 0 at k = {kk:.0e}: A {res['A']:.4f}  B {res['B']:.4f}  LCDM {res['C']:.4f}"
               f"  (A/B - 1 = {res['A']/res['B'] - 1:+.3f})")
txt = "\n".join(out); print(txt); open("iter50_horizon_growth.txt", "w").write(txt + "\n")
