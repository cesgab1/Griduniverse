import numpy as np
G, c, hb, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.3807e-23
lP = np.sqrt(hb*G/c**3); TP = np.sqrt(hb*c**5/G)/kB; rhoP = c**5/(hb*G**2); pP = rhoP*c**2
out = []
# Q1
ls = 1.67*lP
for x in (0.0, 0.5, 0.9):
    P = (lP/ls)**4/(1 - x)
    out.append(f"Q1 grain-set amplitude at T/T_s = {x}: P_Phi ~ {P:.2f}  (measured ~{(3/5)**2*2.1e-9:.1e}) -> x{P/((3/5)**2*2.1e-9):.0e} too large")
# Q2: one free number. Tilts.
ns = 0.9649; nt = 1 - ns
out.append(f"Q2 with one adjusted number: scalar tilt fitted to {ns}; predicted tensor tilt n_t = 1 - n_s = +{nt:.3f} (inflation: n_t = -r/8 <= -0.0045 at r = 0.036)")
# GW background today: Omega_GW(f) ~ 3/128 * Omega_r * r * P_zeta * (f/f*)^n_t  (standard rough transfer for modes entering in radiation era)
Om_r = 9.1e-5; Pz = 2.1e-9; fstar = 0.05/3.086e22*c/(2*np.pi)   # pivot 0.05/Mpc -> Hz
for r in (0.036, 0.001):
    for f, lab in ((1e-3, "LISA"), (25.0, "LIGO")):
        Om = 3/128*Om_r*r*Pz*(f/fstar)**nt
        out.append(f"   r = {r}: Omega_GW at {lab} ({f} Hz) ~ {Om:.1e}  (blue boost x{(f/fstar)**nt:.1f}; LIGO O3 limit ~1.7e-8, LISA reach ~1e-12)")
# Q3: pressure of the hot phase at the fitted temperature (iteration 107: 4.6e-5 T_P) vs grid tensions
g = 106.75; T = 4.6e-5*TP
p = (np.pi**2/90)*g*(kB*T)**4/(hb*c)**3
out.append(f"Q3 pressure to hold at T = {T:.1e} K: {p:.1e} Pa; grid tension today 5.3e-10 Pa (x{p/5.3e-10:.0e} short);"
           f" grid's natural Planck stiffness ~{pP:.1e} Pa (enough). Needed relaxation from then to now: x{pP/5.3e-10:.0e}")
txt = "\n".join(out); print(txt); open("iter108_grid_static.txt", "w").write(txt + "\n")
