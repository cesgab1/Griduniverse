"""
Dark fountain test (follow-up to isolated_well.py): fluid heated on falling into a well, thrown out, must lose its heat
outside (cooling time t_c) so the cosmic fluid stays cold for the Lyman-alpha forest.
Two requirements on t_c:
 (1) Lyman-alpha at z = 3: hot fluid fraction f_hot = (ejection rate per unit fluid) x t_c <~ 1%  (power loss ~ 2 f_hot <~ 2%).
     Ejection rate >= rate at which fresh fluid enters wells deeper than v_th (Press-Schechter, isolated_well.py), ignoring
     recycling (so this is the most lenient t_c,max).
 (2) Escape: the heated fluid must get beyond the galaxy's reach before it cools, or it falls back and the galaxy keeps it
     all (a cycling cloud of size ~ sigma_h t_c holding the full cosmic share -> SPARC/KiDS fail). Lenient reach = turnaround,
     ~3.5 virial radii; MOND's reach is larger (external-field radius V^2/g_e).
     -> sigma_h >= reach / t_c.
 (3) Clusters must keep their fluid: heated_halo.py needs sigma_h <~ 500-600 km/s (fraction kept at 1e13.5-1e15 ~ 1).
"""
import numpy as np
exec(open("isolated_well.py").read().split('print("Fraction')[0])
Gyr = 3.156e16; kpc = 3.086e16                                      # s, km
def age(z): return quad(lambda a: 1/(a*100*h*np.sqrt(Om/a**3+1-Om)), 0, 1/(1+z))[0]*3.086e19/Gyr
def fth(v, z): return erfc(1.686/(np.sqrt(2)*sigM(Mvir(v, z))*D(z)))
z = 3.0; Hz = 100*h*np.sqrt(Om*(1+z)**3+1-Om)/1e3                   # km/s/kpc
ge = 0.025*1.2e-10*3.086e19/1e6                                     # (km/s)^2/kpc
print("Lyman-alpha (z = 3) vs escape, for a typical z = 3 galaxy with circular speed V = 150 km/s")
Rvir = 150/(10*Hz); reach_ta = 3.5*Rvir; reach_mond = 150**2/ge
print(f"  virial radius {Rvir:.0f} kpc; reach: turnaround {reach_ta:.0f} kpc, MOND external-field radius {reach_mond:.0f} kpc")
print("  v_th | ejection rate (/Gyr) | t_c max (Gyr) | sigma_h needed: turnaround / MOND reach (km/s) | clusters allow")
for v in (20, 40, 80, 150):
    rate = (fth(v, 2.8) - fth(v, 3.2))/(age(2.8) - age(3.2))
    tc = 0.01/rate
    s1 = reach_ta*kpc/(tc*Gyr); s2 = reach_mond*kpc/(tc*Gyr)
    print(f"  {v:4d} | {rate:19.3f} | {tc:13.3f} | {s1:8.0f} / {s2:8.0f}{'':24s}| <~ 500-600")
print("\nToday (z = 0), Milky-Way-like galaxy V = 200 km/s:")
H0k = 100*h/1e3; Rv0 = 200/(10*H0k)
for v in (20, 80):
    rate = (fth(v, 0.0) - fth(v, 0.2))/(age(0.0) - age(0.2)); tc = 0.01/rate
    print(f"  v_th {v}: (if the same 1% applied) t_c max {tc:.2f} Gyr -> sigma_h >= {3.5*Rv0*kpc/(tc*Gyr):.0f} km/s (turnaround {3.5*Rv0:.0f} kpc)")
