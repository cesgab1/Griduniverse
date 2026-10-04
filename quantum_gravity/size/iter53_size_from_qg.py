"""
ITERATION 53 (pre-registered in PREREG_53.md): size of dark energy from the QG piece (cells a = 0.49 sqrt(N) l_P, N layers,
spacetime-random sprinkling, one conserved stamp fixed by the whole book's 4-volume tally).
"""
import numpy as np
from scipy.integrate import quad
H0 = 67.4e3/3.0857e22; Om = 0.315; Or = 9.1e-5; OL = 1 - Om - Or; Gyr = 3.156e16; c = 2.998e8
lP, tP, rhoP = 1.616e-35, 5.391e-44, 5.155e96
E = lambda a: np.sqrt(Om/a**3 + Or/a**4 + OL)
Gpc = 3.0857e25
Vvis = 4/3*np.pi*(14.3*Gpc)**3                                   # comoving visible volume today (m^3)
past = quad(lambda a: a**2/(H0*E(a)), 1e-10, 1, limit=200)[0]      # integral a^3 dt (s)
N4 = Vvis*c*past/(lP**3*c*tP)                                     # Planck 4-cells in the visible past
rho_meas = OL*3*H0**2/(8*np.pi*6.674e-11)
r0 = (rhoP/np.sqrt(N4))/rho_meas
out = ["ITERATION 53: size of dark energy from the QG piece (expectations pre-registered)", "",
       f"Planck-cell tally of the visible past N4 = {N4:.2e}; rho_P/sqrt(N4) = {rhoP/np.sqrt(N4):.2e} kg/m^3 = {r0:.2f} x measured "
       f"({rho_meas:.2e})", ""]
f_ind = lambda N: 1/(0.491**2*np.sqrt(N)); f_sh = lambda N: 1/(0.491**2*N)
Smin = (27.5*Gpc)**3/Vvis
out.append(f"CMB topology floor: cube side >= 27.5 Gpc -> S >= {Smin:.2f}")
for cr in (1.0, 2.0):
    Nind = (r0*cr/0.491**2)**2/Smin; Nsh = r0*cr/0.491**2/np.sqrt(Smin)
    out.append(f"c_rule = {cr}: whole book needs S x T = {(r0*cr)**2:.2f} x 17.2/N (independent layers) -> N <= {Nind:.0f};"
               f"  shared grid -> N <= {Nsh:.1f}")
def time_left(T):
    if T < 1: return np.nan
    af = 1.0
    while (past + quad(lambda a: a**2/(H0*E(a)), 1, af)[0])/past < T: af += 0.001
    return quad(lambda a: 1/(a*H0*E(a)), 1, af)[0]/Gyr
out += ["", " c_rule  layers N   S x T   cube side if T=1 [Gpc]   time left if S = floor [Gyr]   cell a [l_P]   cutoff [GeV]"]
EP = 1.22e19
for cr in (1.0, 2.0):
    for N in (4, 10, 20, 44, 100, 118, 176):
        ST = (r0*cr*f_ind(N))**2
        if ST < Smin: out.append(f"  {cr:3.1f}    {N:5d}    {ST:6.2f}   excluded (whole book smaller than the space we already see)"); continue
        L = (ST*Vvis)**(1/3)/Gpc; tl = time_left(ST/Smin); a = 0.491*np.sqrt(N)
        out.append(f"  {cr:3.1f}    {N:5d}    {ST:6.2f}        {L:6.1f}                     {tl:6.1f}                {a:5.2f}       {EP/a:.1e}")
out += ["", "Comparison (not fitted): the Standard Model has ~118 field components (g* = 106.75 counting fermions as 7/8).",
        "GRB check: cells <= 6.5 l_P (= 1e-34 m) vs the measured cell ceiling 5.7e-28 m -> no photon-delay signal expected at any",
        "foreseeable energy (null prediction; a detection would exclude this route)."]
txt = "\n".join(out); print(txt); open("iter53_size_from_qg.txt", "w").write(txt + "\n")
