"""Quark nuggets as ALL of dark matter: what mass range is left? Direct calculations + literature limits marked (mem)."""
import numpy as np
Msun = 1.989e30; mp = 1.673e-27; rho_nuc = 4e17
rho_loc = 0.4*1.783e-27/1e-6          # kg/m^3 local dark matter
v = 2.5e5                              # m/s
R_E = 6.371e6; yr = 3.156e7
print(f"{'nugget mass (kg)':>17s} {'(Msun)':>9s} {'size':>10s} {'quarks/3 (baryon number)':>25s} {'hits on Earth':>22s}")
for M in (1e-3, 1, 1e3, 1e6, 1e9, 1e12, 1e15, 1e18, 1e21, 1e24, 1e27):
    R = (3*M/(4*np.pi*rho_nuc))**(1/3)
    rate = rho_loc*v/M*np.pi*R_E**2*yr
    size = f"{R*1e6:.1f} um" if R < 1e-3 else (f"{R*100:.1f} cm" if R < 10 else f"{R:.0f} m")
    hits = f"{rate:.1e} per yr" if rate >= 1 else f"1 per {1/rate:.1e} yr"
    print(f"{M:17.0e} {M/Msun:9.0e} {size:>10s} {M/mp:25.1e} {hits:>22s}")
print("""
Limits on compact dark lumps as ALL dark matter (mem -- from memory of the literature, not looked up here):
  * light end: lumps hitting Earth often would leave tracks in ancient mica / detectors / seismic data -> excluded up to roughly
    tens of grams to ~1e4 kg (various searches; 'macro' dark matter reviews)
  * microlensing of stars in Andromeda (HSC) and the Milky Way (EROS, OGLE, MACHO): excludes ~2e19 kg (1e-11 Msun) up to ~10 Msun
    for ANY compact object, black hole or nugget
  * open band: roughly 1e4 kg - 2e19 kg (planetoid/asteroid range)
  * survival through the hot early universe (boiling-off of nucleons, Alcock & Farhi 1985; later work with re-absorption): needs
    baryon number above ~1e42-1e52 depending on assumptions -> mass above ~1e15-1e25 kg  [highly uncertain]
Combined: alive only if survival works at baryon number <~ 1e46 (mass <~ 2e19 kg), i.e. within the asteroid band, or if the
'loop' (re-absorption) is efficient. Not excluded, not favoured.""")
