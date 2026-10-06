"""Quantum limits on what the missing pull can be made of (inputs flagged: (v) verified in project, (m) from memory)."""
import numpy as np
h, hbar, c, eV, MS, PC, G = 6.626e-34, 1.0546e-34, 2.998e8, 1.602e-19, 1.989e30, 3.086e16, 6.674e-11
out = []
# systems: name, dark density (Msun/pc^3), speed spread (km/s), size of the dark core (pc)
sy = [("Fornax dwarf (sigma 11.7 (v); core density ~0.1 Msun/pc^3, core ~500 pc (m))", 0.1, 11.7, 500),
      ("Coma at 100 kpc (gas-measured mass, iteration 113 (v))", 7.25e12/(4/3*np.pi*1e5**3), 1000, 1e5)]
out.append("1. If it is made of FERMIONS (like electrons/neutrinos): Pauli exclusion limits how densely they can pack ->")
for n, rho, s, r in sy:
    rho_SI = rho*MS/PC**3; sig = s*1e3
    m4 = rho_SI*h**3/(2*(2*np.pi)**1.5*sig**3)          # max packing: g = 2, occupancy 1
    m = m4**0.25*c**2/eV
    out.append(f"   {n}: particle mass >= {m:.1f} eV")
out.append("2. If it is made of BOSONS (like photons/Higgs): its quantum wave (de Broglie length) must fit inside the core ->")
for n, rho, s, r in sy:
    m = hbar/(s*1e3*r*PC)*c**2/eV
    out.append(f"   {n}: particle mass >= {m:.1e} eV  (lighter -> one quantum wave bigger than the galaxy core)")
out.append("3. It must not stick: Bullet Cluster limit sigma/m < ~1 cm^2/g (m). Ordinary stuff for comparison:")
mp = 1.6726e-24
for n, sgm in [("hydrogen atom on hydrogen atom (~1e-15 cm^2) (m)", 1e-15), ("neutron on proton at ~50 keV (~5 barn) (m)", 5e-24),
               ("neutrino on proton at ~50 keV (~1e-43 cm^2) (m)", 1e-43)]:
    out.append(f"   {n}: {sgm/mp:.1e} cm^2/g -> {'too sticky' if sgm/mp > 1 else 'allowed'}")
out.append("4. If it is black holes (quantum Hawking evaporation): lighter than ~5e11 kg would have evaporated by now (m)")
out.append(f"   -> a black-hole form must weigh >= 5e11 kg (= {5e11*c**2/eV:.0e} eV); the window runs from ~1e-22 eV (bosons) to black holes: >70 orders of magnitude")
txt = "\n".join(out); print(txt); open("portrait.txt", "w").write(txt + "\n")
