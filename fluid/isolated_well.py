"""
User idea: fluid inside a well is cut off from the universe's fluid, so heat gained there stays there; fluid outside stays
cold (and would cool anything mixed into it). Physics translation: heat the fluid LOCALLY when it falls into a well
(fixed heat per unit mass, e.g. latent heat of a state change, triggered by infall faster than v_th), not globally.
Then: shallow wells lose the hot fluid, deep wells keep it (heated_halo.py), and nothing changes with redshift
(lensing at z 0.3 and 0.9 look alike -> passes lensing_z_test by construction), and the large-scale fluid stays cold.
Question tested here: how much of the cosmic fluid passes through wells deeper than v_th by z = 3 (Lyman-alpha) and z = 0?
That is the fluid that gets heated and ejected; unless it cools again outside (the user's 'heat transfer'), it is a hot
component, and the Lyman-alpha forest allows only ~2% small-scale power loss at z = 3.
Press-Schechter mass fraction in halos with circular speed above v_th (BBKS spectrum, sigma8 = 0.81).
"""
import numpy as np
from scipy.special import erfc
from scipy.integrate import quad
Om, h, s8 = 0.31, 0.68, 0.81; rhom = 2.775e11*h**2*Om            # Msun/Mpc^3
def T(k): q = k/(Om*h**2); return np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25
W = lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3
def sig_raw(R): return np.sqrt(quad(lambda lk: (np.exp(lk))**3*np.exp(lk)**0.965*T(np.exp(lk))**2*W(np.exp(lk)*R)**2, -9, 6, limit=400)[0])
norm = s8/sig_raw(8/h)
sigM = lambda M: norm*sig_raw((3*M/(4*np.pi*rhom))**(1/3))
def D(z):
    g = lambda a: quad(lambda x: 1/(x*np.sqrt(Om/x**3+1-Om))**3, 1e-6, a)[0]*2.5*Om*np.sqrt(Om/a**3+1-Om)
    return g(1/(1+z))/g(1.0)
def Mvir(v, z):                                                   # halo mass with circular speed v (km/s), Delta = 200 rho_crit(z)
    H = 100*h*np.sqrt(Om*(1+z)**3+1-Om)/1e3                       # km/s/kpc
    return v**3/(10*4.30e-6*H)
print("Fraction of all matter (hence of the cosmic fluid) in halos with circular speed above v_th:")
print("  v_th (km/s) |   z = 5   z = 3   z = 2   z = 1   z = 0")
for v in (20, 40, 80, 150):
    fr = [erfc(1.686/(np.sqrt(2)*sigM(Mvir(v, z))*D(z))) for z in (5, 3, 2, 1, 0)]
    print(f"  {v:10d}  | " + "  ".join(f"{f:6.2f}" for f in fr))
