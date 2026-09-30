"""
Step 2: what sets a0 when the superfluid is a network of 1D channels?
Superfluid-MOND (Berezhiani-Famaey-Khoury 2018, eqs 3,4,9):  L = (2 Lam (2m)^{3/2}/3) X|X|^{1/2},  L_int = alpha Lam phi rho_b / M_Pl,
a0 = alpha^3 Lam^2 / M_Pl (reduced Planck mass). Fiducial: m = 1 eV, Lam = 0.05 meV, alpha = 5.7 -> a0 = 0.87e-10 m/s^2.
Channel network (Tonks limit): P_1D = (2/3) sqrt(2m) mu^{3/2} / (pi hbar); P_3D = n_ch * P_1D (n_ch = channels crossing unit area)
   => Lam = n_ch / (2 pi m)  (hbar = c = 1).  For spacing l: n_ch ~ 1/l^2 (isotropic random network: ~1/(3 l^2)).
"""
import numpy as np
hbar_eVs, c = 6.582e-16, 2.998e8; eVinv_m = 1.9733e-7; Mp = 2.435e27
acc = lambda E: E/hbar_eVs*c                       # eV -> m/s^2
a0 = lambda al, Lam: al**3*Lam**2/Mp
print(f"check BFK fiducial: a0 = {acc(a0(5.7, 5e-5)):.3e} m/s^2 (paper 0.87e-10)")
H0 = 67.5e3/3.0857e22; target = c*H0/6; tE = target/c*hbar_eVs
rhoL = 3*0.685*(H0*hbar_eVs)**2*Mp**2; lDE = rhoL**-0.25
print(f"dark-energy density^(1/4) = {rhoL**0.25*1e3:.2f} meV -> dark-energy length = {lDE*eVinv_m*1e6:.0f} micron")
for f, lab in ((1.0, "n_ch = 1/l^2"), (1/3, "n_ch = 1/(3 l^2), random isotropic")):
    l = np.sqrt(f/(2*np.pi*1.0*5e-5))
    print(f"BFK fiducial (m=1 eV, Lam=0.05 meV), {lab}: channel spacing l = {l*eVinv_m*1e6:.1f} micron = {l/lDE:.2f} x dark-energy length")
print("\nIf the channel spacing IS the dark-energy length (grid cell = tension scale), what does a0 = cH0/6 demand?")
print("   m (eV)   Lam (meV)   alpha needed")
for m in (0.1, 0.3, 1, 3, 10):
    Lam = 1/(3*2*np.pi*m*lDE**2)
    al = (tE*Mp/Lam**2)**(1/3)
    print(f"   {m:5.1f}   {Lam*1e3:9.4f}   {al:9.1f}")
print("\nUnknowns in a0: particle mass m, channel spacing l, coupling alpha. a0 fixes one combination:")
print("   a0 = alpha^3 n_ch^2 / (4 pi^2 m^2 M_Pl)   -> not derived; needs m and alpha from elsewhere.")
