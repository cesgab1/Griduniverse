"""
Timing check for the 'current' picture: can the fluid stream past forming galaxies/groups yet be trapped by clusters?
Two sources of relative speed between fluid and a forming halo:
 (1) gravity-driven cosmic-web currents: matter (fluid included) falls onto a halo at ~ its own circular speed V_vir;
     anything inside the halo's turnaround radius is bound by construction (spherical collapse), so the current delivers fluid
     INTO halos of every mass (in simulations, filament 'cold streams' penetrate galaxy halos at z~2: Dekel+2009).
 (2) fluid locked to the cosmic rest frame (Khronon's preferred slicing): halos move through it at their peculiar speed
     v_pec(z) (linear theory, rms ~ 300-400 km/s today, lower in the past: v ~ f(z) H(z) a D(z) / (H0 f0 D0) x v0).
     Fluid stays out if v_pec > V_esc of the halo (escape speed at its virial radius).
Need: galaxies AND groups (data: no extra mass) let it pass; clusters (>=1e14) trap it.
"""
import numpy as np
G, Msun, kpc, Mpc = 6.674e-11, 1.989e30, 3.0857e19, 3.0857e22
H0 = 70e3/Mpc; Om = 0.3
Hz = lambda z: H0*np.sqrt(Om*(1+z)**3 + 1 - Om)
def Om_z(z): return Om*(1+z)**3/(Om*(1+z)**3 + 1 - Om)
def D(z):
    a = np.linspace(1e-4, 1/(1+z), 3000); E = np.sqrt(Om/a**3 + 1 - Om)
    return 2.5*Om*E[-1]*np.trapezoid(1/(a*E)**3, a)
f = lambda z: Om_z(z)**0.55
v0 = 350.0                                     # km/s, rms halo peculiar speed today (typical measured/LCDM value)
vpec = lambda z: v0*f(z)*Hz(z)/(1+z)*D(z)/(f(0)*Hz(0)*D(0))
def Vvir(M, z):
    rho = 200*3*Hz(z)**2/(8*np.pi*G); r = (3*M*Msun/(4*np.pi*rho))**(1/3); return np.sqrt(G*M*Msun/r)/1e3
print(f"{'system':28s} {'form z':>6s} {'V_vir':>6s} {'V_esc':>6s} {'v_pec':>6s}  (1) gravity current   (2) frame-locked fluid")
for M, z, lab in ((1e11, 1.1, "small galaxy"), (1e12, 0.9, "Milky-Way-mass galaxy"), (1e13, 0.7, "group"), (3e13, 0.6, "large group"),
                  (1e14, 0.55, "cluster"), (1e15, 0.4, "massive cluster")):
    Vv = Vvir(M, z); Ve = np.sqrt(2)*Vv; vp = vpec(z)
    lock = "passes through" if vp > Ve else "TRAPPED"
    print(f"{lab:28s} {z:6.2f} {Vv:6.0f} {Ve:6.0f} {vp:6.0f}  {'TRAPPED':20s}  {lock}")
print("\nNeeded: galaxies and groups 'passes', clusters 'TRAPPED'.")

# Quantitative: speeds spread as a Maxwellian (1D sigma = v_pec/sqrt 3): trapped fraction = P(relative speed < V_esc).
# Compare with the fraction of the cosmic share (5.4 x baryons) the data need at R500 (gas_retention_test: MOND shortfall / (5.4 x baryon fraction)).
from scipy.special import erf
def trapped(Ve, vrms):
    x = Ve/(vrms/np.sqrt(3)); return erf(x/np.sqrt(2)) - np.sqrt(2/np.pi)*x*np.exp(-x*x/2)
need = {13.0: (0.99, 0.097), 13.5: (0.76, 0.083), 14.0: (0.66, 0.089), 14.5: (0.64, 0.116), 15.0: (0.66, 0.167)}
print(f"\n{'log M500':>8s} {'trapped (frame-locked)':>24s} {'needed by data (fraction of cosmic share)':>44s}")
for lM, (mond, fb) in need.items():
    M = 10**lM; z = 0.7 if lM < 13.8 else 0.5
    Ve = np.sqrt(2)*Vvir(M*1.4, z); tr = trapped(Ve, vpec(z))            # M200 ~ 1.4 M500
    nd = max(1 - mond, 0)/(5.4*fb)
    print(f"{lM:8.1f} {tr:24.2f} {nd:44.2f}")
print("(needed = (1 - MOND/observed) / (5.4 x baryon fraction); hydrostatic-mass bias of ~20% would raise the cluster values)")
