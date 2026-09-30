"""
Test 1: early galaxy formation with grid slack.
Rule (from the Cold Spot analysis): expanding regions feel the cosmic background strain -> no slack; once a region turns around and
falls in (bound), slack switches on: inward pull g_N * nu(g_N/a0).  Spherical collapse of a shell of mass M in a matter-dominated
background (z > 5, dark energy negligible). Find the initial overdensity needed to collapse (reach half the turnaround radius) at z_c,
with and without slack -> ratio of collapse thresholds -> slack threshold = 1.686 x ratio.
Abundance of halos above mass M (Press-Schechter) with sigma(M, z) from CAMB (standard linear growth; slack is off while expanding).
"""
import numpy as np, camb
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import erfc
G = 6.674e-11; Msun = 1.989e30; Mpc = 3.0857e22; a0 = 1.15e-10
h = 0.675; H0 = 100*h*1e3/Mpc; Om = 0.31
rho0 = Om*3*H0**2/(8*np.pi*G)
nu = lambda y: 1/(-np.expm1(-np.sqrt(y)))
OL = 1 - Om
def t_of_a(a): return 2/(3*H0*np.sqrt(OL))*np.arcsinh(np.sqrt(OL/Om)*a**1.5)
def a_of_t(t): return (np.sinh(1.5*H0*np.sqrt(OL)*t)*np.sqrt(Om/OL))**(2/3)
LAM = 3*H0**2*OL
def collapse_a(M, delta_i, a_i, slack):
    # shell starting with Hubble flow (growing mode approx: velocity reduced by delta/3)
    r_i = (3*M/(4*np.pi*rho0*a_i**-3*(1 + delta_i)))**(1/3); t_i = t_of_a(a_i); Hi = H0*np.sqrt(Om*a_i**-3 + OL)
    v_i = Hi*r_i*(1 - delta_i/3)
    state = {"rta": None}
    def f(t, y):
        r, v = y; g = G*M/r**2
        if slack and v < 0: g = g*nu(g/a0)
        return [v, -g + LAM/3*r]
    def ev_ta(t, y): return y[1]
    ev_ta.terminal = False; ev_ta.direction = -1
    sol = solve_ivp(f, [t_i, t_of_a(3.0)], [r_i, v_i], events=ev_ta, rtol=1e-8, atol=1e-6*r_i, dense_output=True, max_step=t_of_a(3.0)/3000)
    if len(sol.t_events[0]) == 0: return None
    t_ta = sol.t_events[0][0]; r_ta = sol.sol(t_ta)[0]
    def ev_c(t, y): return y[0] - r_ta/2
    ev_c.terminal = True; ev_c.direction = -1
    sol2 = solve_ivp(f, [t_ta, t_of_a(3.0)], sol.sol(t_ta), events=ev_c, rtol=1e-8, atol=1e-6*r_i, max_step=(t_of_a(3.0)-t_ta)/4000)
    if len(sol2.t_events[0]) == 0: return None
    tc = sol2.t_events[0][0]; return a_of_t(tc)
def threshold(M, zc, slack, a_i=1e-3):
    ac = 1/(1+zc)
    fn = lambda d: (collapse_a(M, d, a_i, slack) or 10) - ac
    return brentq(fn, 1e-4, 0.5, xtol=1e-7)
zs = [0.0]
pars = camb.set_params(H0=67.5, ombh2=0.0224, omch2=0.12, As=2.1e-9, ns=0.965, mnu=0.06); pars.set_matter_power(redshifts=zs, kmax=300)
res = camb.get_results(pars)
Ms = np.array([1e13, 1e14, 1e15])
Rs = (3*Ms*Msun/(4*np.pi*rho0))**(1/3)/Mpc*h
Rout, zout, sig = res.get_sigmaR(Rs, hubble_units=True, return_R_z=True)
S = {round(z, 1): sig[k] for k, z in enumerate(zout)}
print(f"{'z':>4} {'halo mass':>10} | {'threshold std':>13} {'with slack':>10} | {'sigma':>6} | {'more halos above this mass (x)':>30}")
for zc in (0.0,):
    for j, M in enumerate(Ms):
        dN = threshold(M*Msun, zc, False); dS = threshold(M*Msun, zc, True)
        dc = 1.686*dS/dN; sg = S[round(zc, 1)][j]
        ratio = erfc(dc/np.sqrt(2)/sg)/erfc(1.686/np.sqrt(2)/sg)
        print(f"{zc:4.0f} {M:10.0e} | {1.686:13.3f} {dc:10.3f} | {sg:6.3f} | {ratio:30.2f}")
