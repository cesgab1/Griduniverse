# Where does the electron's lost mass-energy go when it relaxes from +0.78% to today's value?
import numpy as np
h=0.68; Om=0.31; Or=4.15e-5/h**2; OL=1-Om-Or; beta=0.63; ombh2=0.0224; Yp=0.245
rho_crit_eV=1.0537e4*h**2            # eV/cm^3
rho_DE0=OL*rho_crit_eV
n_b0=ombh2*1.878e-29/1.6726e-24      # nucleons /cm^3
n_e0=n_b0*(1-Yp/2)                   # all electrons (H + He)
rho_g0=0.2606                        # CMB photons today, eV/cm^3
def E(z):
    a=1/(1+z); e=np.sqrt(Om*a**-3+Or*a**-4+OL)
    for _ in range(50): e=np.sqrt(Om*a**-3+Or*a**-4+OL*(a*e)**-beta)
    return e
print(f"n_e today {n_e0:.3e} /cm3, dark energy today {rho_DE0:.0f} eV/cm3")
for dm,lab in [(0.0078,'best'),(0.0034,'-1σ'),(0.0122,'+1σ')]:
    q=n_e0*dm*511e3                   # released energy density, comoving-today units (eV/cm3)
    # (A) grid tension = released energy, constant afterwards (pure Λ)
    zA=(rho_DE0/q)**(1/3)-1
    # (B) released energy = ρ_DE(z_s) under the β law ρ_DE ∝ adot^-β
    from scipy.optimize import brentq
    f=lambda z: q*(1+z)**3 - rho_DE0*(E(z)/(1+z))**-beta
    zB=brentq(f,5,5000)
    print(f"Δm/m={dm:.4f} ({lab}): {dm*511:.1f} keV/electron;  switch z (Λ) = {zA:5.0f} [21cm at {1420.4/(1+zA):4.1f} MHz];"
          f"  switch z (β law) = {zB:5.0f} [{1420.4/(1+zB):4.1f} MHz]")
# Alternatives for the energy at z_s=155
zs=155; q=n_e0*0.0078*511e3*(1+zs)**3
print(f"\nat z={zs}: released {q:.0f} eV/cm3; into photons: Δρ/ρ_γ = {q/(rho_g0*(1+zs)**4):.1e} (FIRAS 6e-5) but 4 keV X-rays ionise gas;")
kT=(2/3)*0.0078*511e3/(1+1+ (Yp/4)/(1-Yp))  # shared among electrons+ions, eV
print(f"into gas heat: T ≈ {kT/8.617e-5:.1e} K -> universe fully ionised from z={zs}: optical depth ≈ {0.054*((1+zs)/8.7)**1.5:.1f} (Planck 0.054±0.007)")
# If cracking were exactly at recombination
zr=1090; need=rho_DE0/(n_e0*(1+zr)**3)/511e3
print(f"if the switch were at recombination (z=1090), paying for dark energy needs only Δm/m = {need:.1e}")
for zs in [25,8]:
    print(f"if the switch were at z={zs} (first stars/reionisation): needs Δm/m = {rho_DE0/(n_e0*(1+zs)**3)/511e3:.2f}")
