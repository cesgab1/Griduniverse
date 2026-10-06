"""GRB 080319B (naked-eye burst, z = 0.937): what its timing says about the grid's graininess, vs the record holders."""
import numpy as np
from scipy.integrate import quad
H0 = 67.4/3.086e19; Om = 0.315
E = lambda z: np.sqrt(Om*(1+z)**3 + 1 - Om)
K = lambda z: quad(lambda x: (1+x)/E(x), 0, z)[0]/H0          # seconds (Jacob-Piran lag factor)
EP = 1.22e19                                                   # GeV
dE = 1e-3                                                      # GeV: gamma (~MeV) vs optical (~eV)
for dt in (2.0, 10.0):
    Eqg = dE*K(0.937)/dt
    print(f"080319B: optical and gamma together within {dt:.0f} s over z = 0.937 -> grid/QG scale > {Eqg:.1e} GeV = {Eqg/EP:.0e} Planck energies")
print("record holders already in the project: GRB 090510 (Fermi) linear > 1.2 Planck energies; GRB 221009A (LHAASO) quadratic > 1.2e12 GeV")
print("-> 080319B is ~1e5x weaker as a graininess test (photons only ~MeV); its value is elsewhere (jets, rates, early universe)")
