"""
Shaken-ocean prediction vs galaxy-galaxy lensing across redshift.
Model: before heating (z > 0.65) a galaxy = MOND-boosted baryons + its cosmic share of fluid (unboosted NFW halo, as in
heated_halo.py). After heating (evaporation time R/sigma_h ~ 0.7 Gyr for 300 kpc at 450 km/s; z 0.65 -> 0.37 is ~2.3 Gyr)
the fluid is gone: lensing = MOND-boosted baryons only (which fits KiDS at z ~ 0.2-0.4, chi2 103).
Pre-heating halo mass: Leauthaud+2012 measure Mh/M* ~ 27 at the pivot (constant z = 0.37-0.88); the fluid is (1 - 0.157) of it.
Compare excess surface density DeltaSigma(R) at z ~ 0.37 (after) and z ~ 0.88 (before) at fixed stellar mass.
Data: Leauthaud+2012 (COSMOS, z 0.2-1): pivot Mh/M* constant ~27, low-mass scaling Mh ~ M*^0.46 'does not evolve
significantly'; Hudson+2015 (CFHTLenS, z 0.2-0.8): peak M*/Mh falls 4.5% (z~0.7) -> 3.4% (z~0.3), i.e. halos get relatively
HEAVIER toward today (driven by red galaxies; blue galaxies constant).
"""
import numpy as np
G = 4.30e-6; a0 = 1.2e-10*3.086e19/1e6; rhoc0 = 128.0
nu = lambda y: 1/(1 - np.exp(-np.sqrt(np.maximum(y, 1e-12))))
r = np.logspace(-1, 3.6, 3000)
def nfw_M(Mh, z):
    Ez2 = 0.31*(1+z)**3 + 0.69; R200 = (3*Mh/(4*np.pi*200*rhoc0*Ez2))**(1/3); c = 8.0/(1+z)**0.4
    m = lambda x: np.log(1+x) - x/(1+x); return Mh*m(c*np.minimum(r, 2*R200)/R200)/m(c)
def dsig(M3):
    rho = np.gradient(M3, r)/(4*np.pi*r**2); R = np.logspace(1, 3, 12)
    Sig = np.array([2*np.trapezoid(np.interp(np.sqrt(Rk**2 + l**2), r, rho, right=0), l) for Rk in R
                    for l in [np.logspace(-2, 3.7, 2000)]])
    Mp = np.array([np.trapezoid(2*np.pi*R[:i+1]*Sig[:i+1], R[:i+1]) if i else 0 for i in range(len(R))])
    Rs = np.logspace(-0.5, 3, 400); Ss = np.interp(np.log(Rs), np.log(R), Sig, left=Sig[0])
    Mp = np.array([np.trapezoid(2*np.pi*Rs[Rs <= Rk]*Ss[Rs <= Rk], Rs[Rs <= Rk]) for Rk in R])
    return R, Mp/(np.pi*R**2) - Sig
print("R in kpc; DeltaSigma ratio after/before heating (fixed stellar mass)")
for Ms in (1e10, 4.5e10, 1.5e11):
    Mb = Ms*r**2/(r + 3*(Ms/5e10)**0.3)**2 + Ms*np.minimum(r, 143)/143       # stars + hot gas (Brouwer+21 model)
    mond = nu(G*Mb/r**2/a0)*Mb
    Mh = 27*4.5e10*(Ms/4.5e10)**0.46 if Ms < 4.5e10 else 27*Ms*(Ms/4.5e10)**1.0
    R, after = dsig(mond); _, before = dsig(mond + 0.843*nfw_M(Mh, 0.88))
    print(f"  M* = {Ms:.1e} (Mh before = {Mh:.1e}): " + "  ".join(f"R={Rk:4.0f}: {a/b:.2f}" for Rk, a, b in zip(R, after, before) if Rk in R[[2, 4, 6, 8, 10]]))
print("\nWhat the data allow: inferred halo mass at fixed M* changes by <~0.1 dex (Leauthaud) or rises ~0.1 dex toward today (Hudson);")
print("DeltaSigma at fixed R scales roughly as Mh^0.6-0.8, so an after/before ratio below ~0.75 is excluded at fixed M*.")
