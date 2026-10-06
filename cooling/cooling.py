"""Where does the heat go? Radiation energy lost to expansion vs dark energy gained, in today's visible universe."""
import numpy as np
c = 2.998e8; Mpc = 3.086e22; H0 = 67.4e3/Mpc
rho_c = 3*H0**2/(8*np.pi*6.674e-11); u_c = rho_c*c**2            # J/m^3
V = 4/3*np.pi*(46.5e9*9.461e15)**3                                 # visible universe today (comoving = physical now)
Om_r, Om_L = 9.1e-5, 0.69
E_r, E_L = Om_r*u_c*V, Om_L*u_c*V
yr = 3.156e7
print(f"radiation (CMB + neutrinos) energy in the visible universe now: {E_r:.1e} J;  dark energy: {E_L:.1e} J")
print(f"rates today: radiation LOSES {E_r*H0*yr:.1e} J/yr (each photon stretched); dark energy GAINS {3*E_L*H0*yr:.1e} J/yr (same density, more volume)")
print(f"gain/loss ratio today: {3*Om_L/Om_r:.0f}")
print(f"CMB photons since the fog cleared (z = 1090): energy per photon down x1091 -> {1-1/1091:.2%} of their energy 'gone'")
a_match = E_r/E_L
print(f"cumulative radiation loss equals today's dark energy if counted from a = {a_match:.1e} (z ~ {1/a_match:.0f}); from the grid cap"
      f" (a ~ 2e-32) the loss is {E_r/2e-32:.0e} J -- arbitrary, depends on the start")
print("timing: radiation losses are concentrated EARLY (rate ~ 1/a), dark energy accumulates LATE (~ a^3); CMB limits early dark"
      " energy to ~<1-2% near z ~ 1000 -> 'heat feeds dark energy' fails on timing")
