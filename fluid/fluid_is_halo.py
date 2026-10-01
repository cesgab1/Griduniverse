"""
Test of the one untested idea in Family 4 (IDEAS_LEDGER.md): the dark fluid in a galaxy IS its MOND halo.
The extra mass that MOND implies (the 'phantom' mass, M_dyn - M_b) is real fluid, gathered from the cosmic fluid that the
CMB needs. Not added on top of MOND: one thing. Lensing = dynamics in Khronon, so the lensing mass is that fluid.
Two budget checks, no free numbers except the external field g_e (literature range 0.01-0.05 a0):
 A. Per galaxy, time: phantom mass a galaxy must hold (KiDS lensing follows MOND to ~1 Mpc for isolated lenses, and is the
    same at z ~ 0.9 as at z ~ 0.3), vs the most fluid that can have fallen onto it by then. Spherical shells starting in the
    Hubble flow at z = 200, pulled by the background (matter + Lambda) plus the galaxy's MOND pull (generous: no external
    field, so the pull reaches as far as possible; the fluid's own extra gravity helps only once it is gathered).
 B. Whole universe: total phantom mass of all galaxies (GAMA stellar mass function, Baldry+2012; cold gas and hot gas as in
    our KiDS model) vs the fluid there is (Omega_c = 0.26).
Units: kpc, Msun, km/s, Gyr.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
G = 4.30e-6; a0 = 1.2e-10*3.086e19/1e6; kms_kpc_Gyr = 1.0227      # (km/s)/kpc in 1/Gyr
h = 0.68; H0 = 100*h/1e3*kms_kpc_Gyr; Om, OL = 0.31, 0.69; Ob = 0.049; Oc = Om - Ob
rhoc = 3*(H0/kms_kpc_Gyr)**2/(8*np.pi*G)                          # Msun/kpc^3
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def tz(z): return quad(lambda a: 1/(a*H0*np.sqrt(Om/a**3 + OL)), 0, 1/(1+z))[0]
G_ = G*kms_kpc_Gyr**2                                              # kpc^3/(Msun Gyr^2)
a0_ = a0*kms_kpc_Gyr**2
def collapsed_by(R_L, Mb, t_end, ai=1/201):
    """Shell at comoving radius R_L (kpc): does it fall to r < 1% of start-of-collapse scale by t_end?"""
    r0 = ai*R_L; Hi = H0*np.sqrt(Om/ai**3 + OL)
    # background acceleration of a comoving shell: -(4 pi G/3) rho_m(t) r + (Lambda) r; do it with the scale factor
    def ev(t, y): return y[0] - 0.01*r0
    ev.terminal = True
    def rhs(t, y):
        r, v = y
        # mean matter inside the shell is conserved (shell-crossing ignored): M_m = 4/3 pi R_L^3 Om rhoc
        Mm = 4/3*np.pi*R_L**3*Om*rhoc
        gN = G_*Mb/r**2
        g = G_*Mm/r**2 - OL*H0**2*r + nu(gN/a0_)*gN - gN       # Newtonian background + Lambda + MOND extra from the galaxy
        return [v, -g]
    ti = tz(1/ai - 1)
    s = solve_ivp(rhs, [ti, t_end], [r0, Hi*r0], events=ev, rtol=1e-7, atol=1e-9, max_step=0.05)
    return s.status == 1
def reach(Mb, z):
    t = tz(z); lo, hi = 10.0, 2e4
    for _ in range(40):
        mid = np.sqrt(lo*hi)
        if collapsed_by(mid, Mb, t): lo = mid
        else: hi = mid
    return lo
print("A. Per galaxy: fluid that can have fallen in vs phantom mass the lensing needs")
print("   Mb (Msun) | v_flat | gathered by z=0.9 / z=0 | phantom needed within 300 kpc / 1 Mpc | ratio (gathered z=0.9 / needed 300 kpc)")
for Mb in (1e10, 1e11, 3e11):
    vf = (G*Mb*a0)**0.25
    gath = [Oc*rhoc*4/3*np.pi*reach(Mb, z)**3 for z in (0.9, 0.0)]
    need = [vf**2*R/G - Mb for R in (300, 1000)]
    print(f"   {Mb:9.0e} | {vf:5.0f}  | {gath[0]:9.2e} / {gath[1]:9.2e} | {need[0]:9.2e} / {need[1]:9.2e} | {gath[0]/need[0]:.2f}")
print("   (cosmic share for comparison: 5.4 x Mb)")

print("\nB. Whole universe: phantom mass of all galaxies out to the external-field radius, M_ph = Mb (a0/g_e - 1)")
# GAMA double Schechter (Baldry+2012, h = 0.7): M* = 10^10.66, phi1 = 3.96e-3, a1 = -0.35, phi2 = 0.79e-3, a2 = -1.47 per Mpc^3 per dex
lm = np.linspace(7, 12.5, 400); m = 10**(lm - 10.66)
phi = np.log(10)*np.exp(-m)*(3.96e-3*m**(-0.35+1) + 0.79e-3*m**(-1.47+1))*(h/0.7)**3   # per Mpc^3 per dex
Ms = 10**lm; fgas = np.minimum(10**(-0.5*(lm - 9.5)), 10.0)          # cold gas / stars (rises to dwarfs)
Mb = Ms*(1 + fgas) + Ms                                            # + hot CGM = M* (our KiDS model)
rho_b_gal = np.trapezoid(phi*Mb, lm)/1e9                           # Msun/kpc^3
print(f"   baryons in galaxies (stars + cold + hot gas): Omega = {rho_b_gal/rhoc:.4f}  ({rho_b_gal/rhoc/Ob*100:.0f}% of all baryons)")
for ge in (0.01, 0.025, 0.05):
    Om_ph = rho_b_gal*(1/ge - 1)/rhoc
    print(f"   g_e = {ge:5.3f} a0: Omega_phantom = {Om_ph:.3f}  vs fluid available Omega_c = {Oc:.3f}  -> ratio {Om_ph/Oc:.2f}")

# Groups and clusters (>= 1e13 Msun) hold more baryons than galaxies; their gas also sits in deep-MOND outskirts.
# Rough: matter fraction in halos > 1e13 at z = 0 ~ 0.2 (Press-Schechter), baryon fraction inside ~ 0.1 of the halo mass.
Ob_grp = 0.2*Om*0.10
for ge in (0.025, 0.05):
    print(f"   + groups/clusters (Omega_b ~ {Ob_grp:.4f}), g_e = {ge} a0: total ratio ~ {(rho_b_gal/rhoc + Ob_grp)*(1/ge - 1)/Oc:.1f}")
