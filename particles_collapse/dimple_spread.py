"""
Does the electron cloud (or anything physical in the atom) spread each atom's dimple over >= 0.54 Angstrom, the size strain-collapse needs?
The dimple follows MASS. Where is an atom's mass, and how spread out is it?
  - electron cloud: spread ~1 Angstrom, but carries only Z*m_e / M_atom of the mass
  - nucleus: 99.97% of the mass, spread only by its zero-point/thermal jiggle in the crystal (Debye model)
Effective spread for strain collapse = the mass-weighted spread seen in the self-energy (dominated by the nucleus).
Detector: germanium at ~90 K (Gran Sasso), Debye temperature 374 K.
"""
import numpy as np
hbar = 1.0546e-34; kB = 1.3807e-23; u = 1.6605e-27; me = 9.109e-31
M = 72.63*u; Z = 32; thD = 374.0
def urms(T):
    # Debye model, mean-square displacement along one axis
    x = np.linspace(1e-6, thD/T if T > 0 else 1e3, 20000)
    integral = np.trapezoid(x/(np.expm1(x)), x) if T > 0 else 0
    return np.sqrt(3*hbar**2/(M*kB*thD)*(0.25 + (T/thD)**2*integral))
fe = Z*me/M
print(f"germanium: electron cloud carries {fe*100:.3f}% of the atom's mass; nucleus {100-fe*100:.3f}%")
for T in (0, 90, 300):
    print(f"   nucleus jiggle at {T:3d} K: {urms(T)*1e10:.3f} Angstrom (needed: >= 0.54)")
# strain self-energy with a 2-part mass distribution: nucleus Gaussian (width u) + electron cloud Gaussian (width 1 A)
def self_energy(sig_n, sig_e=1e-10):
    # gravitational self-energy of Gaussian blobs ~ G m^2 /(2 sqrt(pi) sigma); cross term with combined width
    fn, fe_ = 1 - fe, fe
    s = lambda a, b: 1/np.sqrt(np.pi*(a**2 + b**2)/2)/2
    return fn**2*s(sig_n, sig_n) + fe_**2*s(sig_e, sig_e) + 2*fn*fe_*s(sig_n, sig_e)
R_eff = 1/(2*np.sqrt(np.pi))/self_energy(urms(90))      # equivalent single-Gaussian width
print(f"\neffective dimple spread of a germanium atom at 90 K (nucleus + electron cloud): {R_eff*1e10:.3f} Angstrom")
print(f"-> {0.54e-10/R_eff:.0f}x smaller than the Gran Sasso minimum: collapse would be ~{(0.54e-10/R_eff):.0f}x faster than allowed")
print("   electron cloud's share of the self-energy:", f"{(fe**2/(1e-10) ) / (self_energy(urms(90))*2*np.sqrt(np.pi)) *100:.4f}%")
