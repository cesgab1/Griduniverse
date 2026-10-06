"""How transparent would a universe full of dense lumps (e.g. quark nuggets) be? Optical depth = column mass x cross-section per mass."""
import numpy as np
rho_nuc = 4e17                       # kg/m^3, nuclear density
rho_dm = 0.26*8.5e-27                 # kg/m^3, average dark matter today
c, H0 = 2.998e8, 2.18e-18
col_universe = rho_dm*c/H0            # kg/m^2 across a Hubble length
col_halo = 0.4*1.783e-27/1e-6*3.086e20*10   # local density 0.4 GeV/cm^3 over ~10 kpc, kg/m^2
for M in (1e14, 1e17, 1e20):          # kg (1e17, 1e20, 1e23 g)
    R = (3*M/(4*np.pi*rho_nuc))**(1/3); s = np.pi*R**2/M
    print(f"lump {M:.0e} kg: radius {R*100:.1f} cm; blocking area per kg {s:.1e} m^2/kg; "
          f"fraction of light blocked across the whole universe {col_universe*s:.0e}, across our galaxy's halo {col_halo*s:.0e}")
