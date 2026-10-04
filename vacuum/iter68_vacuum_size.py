"""ITERATION 68 (pre-registered in PREREG_68.md)."""
import numpy as np
from scipy.integrate import quad
H0 = 67.4e3/3.0857e22; Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or; c = 2.998e8
hbar, G = 1.054571817e-34, 6.674e-11
E = lambda a: np.sqrt(Om/a**3 + Or/a**4 + OL)
t_of = lambda a: quad(lambda x: 1/(x*H0*E(x)), 1e-12, a, limit=300)[0]
chi_of = lambda a: quad(lambda x: c/(x*x*H0*E(x)), a, 1, limit=300)[0]                     # comoving distance from a to today
Vcyl = 4/3*np.pi*(14.3*3.0857e25)**3*c*quad(lambda a: a**2/(H0*E(a)), 1e-12, 1, limit=300)[0]
avals = np.geomspace(1e-8, 1, 3000)
integrand = [x**3*4/3*np.pi*chi_of(x)**3/(x*H0*E(x)) for x in avals]                       # dt = da/(a H)
Vlc = c*np.trapezoid(np.array(integrand)*avals, np.log(avals))
rho_meas = OL*3*H0**2/(8*np.pi*G)                                                          # kg/m^3
out = ["ITERATION 68: vacuum-energy route to the size -- every convention (expectations pre-registered)", "",
       f"4-volumes: past light cone {Vlc:.2e} m^4, visible-region cylinder {Vcyl:.2e} m^4 (ratio {Vcyl/Vlc:.1f})", "",
       " cell length            4-volume         factor   rho_DE / measured"]
lP2 = hbar*G/c**3
res = []
for lname, l2 in (("Planck", lP2), ("reduced Planck", 8*np.pi*lP2)):
    for vname, V in (("past light cone", Vlc), ("visible cylinder", Vcyl)):
        for f in (0.5, 1.0):
            rho = f*hbar*c/(l2*np.sqrt(V))/c**2                                             # kg/m^3
            res.append(rho/rho_meas)
            out.append(f" {lname:15s}   {vname:17s}   {f:4.1f}     {rho/rho_meas:8.3f}")
out += ["", f"spread: {min(res):.3f} - {max(res):.2f} x measured ({np.log10(max(res)/min(res)):.1f} orders of magnitude)",
        "Each value would TRACK the total density if re-set continuously (excluded); it can only be today's size if frozen."]
txt = "\n".join(out); print(txt); open("iter68_vacuum_size.txt", "w").write(txt + "\n")
