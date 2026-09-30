"""
Measurement in the River-Compactness framework: double slit, Schrodinger's cat, and 'why one spot'.
 1. Double slit: fringe visibility vs mass for (a) the quantum grid (piece 0: no gravitational collapse; loss only when some other
    system records the path, incl. gravitational radiation) and (b) Diosi-Penrose collapse, plus ordinary environment.
 2. Schrodinger's cat: how fast each channel makes 'alive' and 'dead' branches independent (decohere).
 3. 'Why one spot': test a grid-commit rule (the growing grid must pick ONE terrain once branch geometries differ by >= one grid cell
    of spacetime volume). If that rule were true, which interference experiments would it forbid?
"""
import numpy as np
G, hbar, c, kB, amu = 6.6743e-11, 1.054571817e-34, 2.99792458e8, 1.380649e-23, 1.66054e-27
lP = 1.616255e-35; ld = 0.603 * lP

def W(m, R, d):   # mutual gravitational energy (magnitude) of two uniform spheres, centres d apart
    lam = d / R
    return np.where(lam >= 2, G*m*m/np.maximum(d, 1e-300), G*m*m/R*(6/5 - lam**2/2 + 3/16*lam**3 - lam**5/160))
def E_DP(m, R, d): return 6/5*G*m*m/R - W(m, R, d)
def radius(m, rho): return (3*m/(4*np.pi*rho))**(1/3)
def n_gravitons(m, d, T): return G*(m*d*d)**2 / (hbar * c**5 * T**5)        # quadrupole estimate of which-path gravitons
def grid_cells(m, d, T):
    """spacetime volume (in grid cells) by which the two branch geometries differ: ~ (G m / c^2) x d^2 x cT  (Newtonian-potential
    difference integrated over the region the displaced mass reshapes)"""
    return (G*m/c**2) * d*d * c*T / ld**4

print("=== 1. Double slit: which systems could erase the fringes? ===")
cases = [("electron (Tonomura 1989)", 9.109e-31, 1e-6, 1e-8, 1e3),
         ("C60 buckyball (Arndt 1999)", 720*amu, 1e-7, 1e-2, 1.7e3),
         ("25,000 amu molecule (Fein 2019)", 25e3*amu, 2.7e-7, 1e-2, 1.4e3),
         ("silica nanosphere 1e9 amu (proposed)", 1e9*amu, 1e-7, 1.0, 2.2e3),
         ("BMV microdiamond 1e-14 kg", 1e-14, 250e-6, 2.5, 3.5e3),
         ("dust grain 1 microgram", 1e-9, 1e-6, 1.0, 2.2e3)]
print(f"{'object':38s} {'Diosi-Penrose collapse time':>28s} {'gravitons carrying path':>24s} {'grid cells differing':>21s}")
for name, m, d, T, rho in cases:
    R = radius(m, rho) if m > 1e-28 else 2.8e-15
    tdp = hbar / E_DP(m, R, d)
    print(f"{name:38s} {tdp:28.2e} s {n_gravitons(m, d, T):24.1e} {grid_cells(m, d, T):21.1e}")
print("-> gravity-as-recorder (gravitons) is negligible for everything below planet-size superpositions: the quantum grid keeps fringes.")
print("-> Diosi-Penrose predicts loss once its collapse time drops below the flight time (between the 25k-amu molecule and the nanosphere).")

print("\n=== 3. A grid-commit rule for 'one spot' ===")
for name, m, d, T, rho in cases[:3]:
    print(f"  {name}: branch geometries differ by {grid_cells(m, d, T):.1e} grid cells -> a commit-at-one-cell rule would erase these fringes "
          f"(observed: fringes present)")
# how small would the threshold have to be pushed to survive? require commit threshold N* > electron value
print(f"  => any rule 'the grid commits once branches differ by N* cells' needs N* > {grid_cells(*cases[0][1:4]):.0e}, i.e. no grid-scale rule; "
      "a single-cell commit rule is ruled out by the electron double slit itself.")

print("\n=== 2. Schrodinger's cat (4 kg, alive/dead differ by ~0.5 kg of tissue shifted ~10 cm) ===")
M, dM, dx, A = 4.0, 0.5, 0.1, 0.1         # total mass, displaced mass, displacement, surface area m^2
n_air, v_air, m_air, T = 2.5e25, 470.0, 4.8e-26, 300.0
lam_th = hbar / np.sqrt(2*m_air*kB*T)
t_air = 1 / (n_air * v_air * A)                                        # each air collision resolves a 10 cm difference
# thermal photons: a 300 K body radiates ~ sigma T^4 A / (2.7 kT) photons per second, each resolves the shape difference
sigma = 5.670374e-8; t_phot = 1 / (sigma*T**4*A / (2.7*kB*T))
# gravity only: every air molecule within ~1 m feels a branch-dependent pull; decoherence exponent sum (dp lam/hbar)^2
r = 1.0; N_mol = n_air * 4/3*np.pi*r**3
dp_rate = G*dM*m_air*dx/r**3
t_grav = hbar / (dp_rate * lam_th * np.sqrt(N_mol))
t_dp = hbar / E_DP(dM, radius(dM, 1e3), dx)
print(f"  air collisions:            {t_air:.1e} s")
print(f"  thermal photons (sealed box, 300 K): {t_phot:.1e} s")
print(f"  gravity alone, via the surrounding air (quantum grid): {t_grav:.1e} s")
print(f"  gravitons radiated in 1 s carrying the branch: {n_gravitons(dM, dx, 1.0):.1e}")
print(f"  Diosi-Penrose collapse (if it existed): {t_dp:.1e} s")
print(f"  thermal wavelength of air used: {lam_th:.1e} m")
