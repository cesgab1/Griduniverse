"""Work backwards from what clusters need (Coalesce): size of the extra pull, and what could supply it if it came from cooling."""
import numpy as np
Ms, c, kB, mp = 1.989e30, 2.998e8, 1.381e-23, 1.673e-27
M_vis, M_tot = 3e14, 1.6e15            # Coma (iteration 94/98 numbers)
M_extra = M_tot - M_vis; E_need = M_extra*Ms*c**2
print(f"Coma: extra pull equivalent to {M_extra:.1e} Msun = {E_need:.1e} J")
# 1 heat in Coma's gas now
N = 2.5e14*Ms/(0.6*mp); E_heat = 1.5*N*8.2e3*1.602e-19
print(f"  heat in Coma's hot gas now: {E_heat:.1e} J -> short x{E_need/E_heat:.0e}")
# 2 all X-ray cooling over cosmic history
L = 8e37; E_cool = L*13.8e9*3.156e7
print(f"  X-ray cooling at today's rate for 13.8 Gyr: {E_cool:.1e} J -> short x{E_need/E_cool:.0e}")
# 3 expansion cooling of the radiation in the region that later fell into Coma
rho_m0 = 0.31*2.775e11*0.674**2                     # Msun/Mpc^3
V = M_tot/rho_m0*(3.086e22)**3                       # m^3, comoving
u_r0 = 9.1e-5*7.65e-10                               # J/m^3 radiation today
for z in (1090, 3400):
    print(f"  expansion cooling of that region's radiation since z = {z}: {u_r0*V*z:.1e} J -> ratio {u_r0*V*z/E_need:.2f}"
          f"  (matches in size because matter-radiation equality is DEFINED by the dark matter amount: circular)")
# 4 shape: extra pull law inferred from data, galaxies vs clusters
a_gal, a_cl = 1.2e-10, 2.0e-9
Tvir = lambda v: 0.6*mp*(v*1e3)**2/(2*kB)
print(f"inferred law in both: extra pull ~ sqrt(g_bar x scale); scale galaxies {a_gal:.1e}, clusters {a_cl:.1e} (x{a_cl/a_gal:.0f})")
print(f"  'temperature' of the systems: galaxy (150 km/s) {Tvir(150):.0e} K, cluster (1000 km/s) {Tvir(1000):.0e} K (x{Tvir(1000)/Tvir(150):.0f})"
      f" -> the HOTTER system needs the LARGER scale: opposite to 'cooling makes the pull'")
