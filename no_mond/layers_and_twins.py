"""
Revisit (Coalesce, Oct 2026): two earlier pictures, judged WITHOUT MOND.
 (1) 'Electron twins' on the mosaic: (a) the fake copies (doublers) a grid gives every electron; (b) a light electron-like
     twin fermion as the dark component.
 (2) 'Fluid / gas between the grid layers': the dark fluid lives in narrow channels between layers; its exact 1-D gas law
     (fluid/channel_eos.txt) is P ~ rho^2 when dense and P ~ rho^3 when dilute. Earlier it was judged as a MOND superfluid.
     Without MOND it is simply a dark fluid with ONE global pressure law. That fixes how a halo's core size depends on its
     central density:
        Gamma = 2 (n = 1 polytrope):   core radius fixed, the same in every galaxy       -> slope d ln r_core / d ln rho_core = 0
        Gamma = 3 (n = 1/2 polytrope): r_core^2 ~ K rho_core                              -> slope +1/2
     Test: fit cored (Burkert) halos to SPARC in Newtonian gravity (no MOND), measure the slope, compare.
     Also: does a pressure law strong enough to make kpc cores spoil the early universe (Jeans length vs Lyman-alpha)?
"""
import numpy as np, glob
from scipy.optimize import minimize
rng = np.random.default_rng(1)
G = 4.30e-6; rhoc = 136.0                                   # kpc (km/s)^2 / Msun ; Msun/kpc^3
fB = lambda x: np.log(1 + x**2) + 2*np.log(1 + x) - 2*np.arctan(x)
def bur_v2(r, M, rc):
    R200 = (3*M/(4*np.pi*200*rhoc))**(1/3); return G*M*fB(r/rc)/fB(R200/rc)/r
def chi(v2, V, e): return np.sum(((np.sqrt(np.maximum(v2, 0)) - V)/e)**2)
rows = []
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    ok = (r > 0) & (V > 0); r, V, eV, Vg, Vd, Vb = r[ok], V[ok], eV[ok], Vg[ok], Vd[ok], Vb[ok]
    if len(r) < 8: continue
    vb2 = Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2; e = np.sqrt(eV**2 + (0.05*V)**2)
    best = None
    for rc0 in (0.5, 2, 8):
        s = minimize(lambda p: chi(vb2 + bur_v2(r, 10**p[0], 10**p[1]), V, e), [11, np.log10(rc0)], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-4, maxiter=800))
        if best is None or s.fun < best.fun: best = s
    M, rc = 10**best.x[0], 10**best.x[1]
    R200 = (3*M/(4*np.pi*200*rhoc))**(1/3); rho0 = M/(np.pi*rc**3*fB(R200/rc))
    # keep only well-measured cores: core inside the data range, decent fit, not a pure-baryon galaxy
    good = (r[1] < rc < 0.8*r[-1]) and best.fun/len(r) < 3 and 1e5 < rho0 < 1e10
    rows.append((fn.split("/")[-1][:-11], rc, rho0, V[-1], good))
names, rcs, rhos, vfl, good = zip(*rows); rcs, rhos, vfl, good = map(np.array, (rcs, rhos, vfl, good))
x, y = np.log10(rhos[good]), np.log10(rcs[good])
sl = np.polyfit(x, y, 1)[0]; bs = []
for _ in range(2000):
    i = rng.integers(0, len(x), len(x)); bs.append(np.polyfit(x[i], y[i], 1)[0])
print("(2) FLUID BETWEEN THE LAYERS, no MOND: SPARC cored-halo fits (Newtonian)")
print(f"    {good.sum()} galaxies with a well-measured core (of {len(rows)})")
print(f"    core radius {np.median(rcs[good]):.1f} kpc median, spread x{10**np.std(y):.1f} (0.5-95%: {np.percentile(rcs[good],5):.1f}-{np.percentile(rcs[good],95):.1f} kpc)")
print(f"    slope d ln r_core / d ln rho_core = {sl:+.2f} +/- {np.std(bs):.2f}   (Donato+2009 surface-density law: -1)")
print(f"    channel gas law predicts:  Gamma=2 -> 0  ({abs(sl)/np.std(bs):.0f} sigma off);  Gamma=3 -> +0.5  ({abs(sl-0.5)/np.std(bs):.0f} sigma off)")
S = rhos[good]*rcs[good]/1e6
print(f"    rho_core x r_core = {np.median(S):.0f} Msun/pc^2 (spread x{10**np.std(np.log10(S)):.1f}) -- the observed near-constant combination")
print(f"    rho_core x r_core^2 (Gamma=3 constant would be rho/r^2) spread x{10**np.std(np.log10(rhos[good]/rcs[good]**2)):.1f};"
      f" fixed r_core (Gamma=2) spread x{10**np.std(y):.1f}")

# early-universe check: pressure strong enough for a 3 kpc core at the median core density
kpc = 3.086e19; Msun = 1.989e30; Gsi = 6.674e-11; c = 2.998e8
rho_core = np.median(rhos[good])*Msun/kpc**3; Rc = np.median(rcs[good])*kpc
rho_bar0 = 0.26*3*(67.7e3/3.086e22)**2/(8*np.pi*Gsi)
print("\n    Early universe (pressure law fixed by the median core; Jeans length must stay < ~0.1 Mpc comoving for Lyman-alpha):")
for Gam, xi1, nn in ((2, np.pi, 1.0), (3, 2.7528, 0.5)):
    K = (Rc/xi1)**2*4*np.pi*Gsi/((nn + 1)*rho_core**(1/nn - 1))   # Lane-Emden radius R = xi1 sqrt((n+1) K rho_c^(1/n-1) / 4 pi G)
    for z in (0, 10, 100, 1100):
        rho = rho_bar0*(1 + z)**3; cs2 = Gam*K*rho**(Gam - 1)
        lamJ = np.sqrt(cs2)*np.sqrt(np.pi/(Gsi*rho/0.84))*(1 + z)/3.086e22   # comoving Mpc, total matter ~ rho_dm/0.84
        print(f"      Gamma={Gam}, z={z:5d}: sound speed {np.sqrt(cs2)/1e3:10.3g} km/s ({np.sqrt(cs2)/c:.1e} c), comoving Jeans length {lamJ:9.3g} Mpc")

# (1b) a fermion 'electron twin' as ALL the dark matter: Pauli (Tremaine-Gunn) limit from the same SPARC cores
h = 6.626e-34; eV = 1.783e-36
sig = vfl[good]/np.sqrt(2)*1e3
m_min = (rhos[good]*Msun/kpc**3*h**3/(2*(2*np.pi*sig**2)**1.5))**0.25/eV
print("\n(1b) ELECTRON-LIKE TWIN FERMION as all the dark matter (no MOND): Pauli limit from SPARC cores")
print(f"    minimum mass: median {np.median(m_min):.0f} eV, max {m_min.max():.0f} eV (densest dwarf cores)")
print("    vs the 2 eV twin of the earlier cluster test: 2 eV cannot hold galaxy halos (Pauli), and CMB allows omega_f < ~0.002 at 2 eV")
print("    heavier twins: Lyman-alpha warm-DM bound m > ~5.7 keV (thermal relic, Irsic+2024) -> then it behaves as cold DM")
