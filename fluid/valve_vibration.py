"""
User ideas (brainstorm): (A) a Tesla-valve-like grid (flow passes easily one way, is blocked the other way);
(B) grid vibrations changing how the fluid behaves (here: shaking/heating the fluid late in cosmic history).

A. Tesla valve, physics facts used: a Tesla valve only rectifies at high Reynolds number (inertial flow); at slow flow it is
   symmetric. So it is a speed switch that ADDS trapping (fluid falls in fast, can't bounce back out). It cannot stop slow,
   ordinary gravitational capture by galaxies, which is the core problem (clusters/README.md lesson). Qualitative only.
   Pumping: a valve plus oscillation = a valveless pump. That makes A a way to turn B's vibrations into flow (kept for later).

B. Late heating by grid vibrations ('shaken ocean'):
   1. Energy source and timing from our dark-energy law rho_DE ~ adot^-1/2: rho_DE rises while the universe decelerates and
      starts FALLING when acceleration begins (adot minimum). Released energy could shake the grid. Timing is derived, not picked.
   2. Requirement: heat the fluid to a speed spread sigma_h so galaxies/groups (V_esc < ~3 sigma_h) lose their halos, clusters keep them.
   3. Shape fit: retained fraction = P(v < V_esc) for a Maxwellian of 1-D spread sigma_h (one number).
   4. Cost: a fluid with spread sigma_h stops clustering below its Jeans length -> total matter growth (sigma8, S8) drops.
"""
import numpy as np
from scipy.special import erf
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import solve_ivp
exec(open("current_timing.py").read().split("print(f\"{'system'")[0])
c = 2.998e5
Om_, Ob_ = 0.31, 0.049; Of = Om_ - Ob_; OL = 1 - Om_

# --- B1: energy released by the dark-energy law and its timing
E = lambda a: np.sqrt(Om_/a**3 + OL)                      # background ~LCDM (our fits are within ~1% of it)
adot = lambda a: a*E(a)
a_min = minimize_scalar(adot, bounds=(0.3, 1), method="bounded").x
z_acc = 1/a_min - 1
rel = 1 - (adot(a_min)/adot(1))**0.5                       # fractional fall of rho_DE since acceleration began
released = rel*OL                                          # in units of rho_crit c^2
print(f"B1. Acceleration starts at z = {z_acc:.2f}; the DE law then releases {rel*100:.1f}% of rho_DE = {released:.3f} rho_crit c^2 by today")

# --- B2/B3: shape fit with the required heating
target = {12.0: 0.0, 13.0: 0.02, 13.5: 0.54, 14.0: 0.71, 14.5: 0.57, 15.0: 0.38}
target_b = {12.0: 0.0, 13.0: 0.50, 13.5: 1.0, 14.0: 1.0, 14.5: 0.97, 15.0: 0.65}   # 20% hydrostatic bias, capped at 1
def kept(Ve, s): x = Ve/s; return erf(x/np.sqrt(2)) - np.sqrt(2/np.pi)*x*np.exp(-x*x/2)
Ve0 = {lM: np.sqrt(2)*Vvir(10**lM*1.4, 0.0) for lM in target}
def fit(tg):
    cost = lambda s: sum(((kept(Ve0[lM], s) - t)/0.15)**2 for lM, t in tg.items())
    r = minimize_scalar(cost, bounds=(50, 3000), method="bounded"); return r.x, r.fun
print("    escape speeds today (km/s): " + "  ".join(f"1e{lM}: {Ve0[lM]:.0f}" for lM in target))
for name, tg in [("as measured", target), ("20% hydrostatic bias", target_b)]:
    s, chi = fit(tg)
    print(f"B3. {name}: sigma_h = {s:.0f} km/s, mismatch chi2 = {chi:.1f} | kept: " +
          "  ".join(f"{kept(Ve0[lM], s):.2f}(need {tg[lM]:.2f})" for lM in tg))
    need = Of*1.5*(s/c)**2
    print(f"    heating energy {need:.1e} rho_crit c^2 = {need/released:.0e} of the DE release (gravitational-wave background ~1e-9: {need/1e-9:.0e}x too small)")

# --- B4: growth cost. Two fluids (baryons cold, dark fluid heated to sigma_h from z_acc on), linear, Newtonian.
def growth(k_hMpc, s):
    k = k_hMpc*0.68*(c/70.0)                                 # k in units of H0/c
    def rhs(lna, y):
        a = np.exp(lna); H = E(a); dlnH = -1.5*Om_/a**3/H**2
        df, dfp, db, dbp = y
        cs2 = (s/c)**2 if a > a_min else 0.0
        src = 1.5/(a**3*H**2)*(Of*df + Ob_*db)
        return [dfp, -(2 + dlnH)*dfp + src - cs2*k**2/(a**2*H**2)*df, dbp, -(2 + dlnH)*dbp + src]
    a0 = 1e-2; y0 = [a0, a0, a0, a0]
    sol = solve_ivp(rhs, [np.log(a0), 0], y0, rtol=1e-8, atol=1e-12)
    df, _, db, _ = sol.y[:, -1]; return (Of*df + Ob_*db)/Om_
# sigma8 with a BBKS spectrum (n_s = 0.965, Gamma = Om h), top-hat R = 8 Mpc/h
ks = np.logspace(-3, 1.3, 120); h = 0.68
q = ks/(Om_*h); T = np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25
P = ks**0.965*T**2; W = lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3
g0 = np.array([growth(k, 0) for k in ks])
def s8(s): g = np.array([growth(k, s) for k in ks]); return np.sqrt(np.trapezoid(P*(g*W(8*ks))**2*ks**2, ks)/np.trapezoid(P*(g0*W(8*ks))**2*ks**2, ks))
for s in (300, 600, 800, 1100, 1500):
    print(f"B4. sigma_h = {s:4d} km/s from z = {z_acc:.2f}: sigma8 / unheated = {s8(s):.3f}")
print("    (weak-lensing S8 sits ~5-8% below Planck: a drop of up to ~8% is welcome, > ~12% excluded)")
print("    Caveat: fluid near halo centres is bound ~2.5x deeper than at R500; if it must leave too, sigma_h rises ~2.5x.")
