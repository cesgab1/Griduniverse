"""Step 1-2: light fermions ('lighter electrons', i.e. sterile-neutrino-like) collecting in clusters.
Need: extra mass f x baryons so that slack acting on (baryons + fermions) gives the measured cluster pull:
      (g_b(1+f)) ν(g_b(1+f)/a0) = g_meas(g_b),  a0 = cH0/6, CLASH/BCG relation (Tian+2024).
Clusters cannot hold more fermions than the cosmic share: ω_f >= f ω_b (fermions collected from the same region as the gas).
Pauli (Tremaine-Gunn) capacity: max density of fermions of mass m with speed spread σ (2 spin states): ρ_max = 2 m^4 (2πσ²)^{3/2}/h^3 (non-degenerate bound)."""
import numpy as np
from scipy.optimize import brentq
c=2.998e8; Mpc=3.0857e22; H0=68.5e3/Mpc; a0=c*H0/6; G=6.674e-11; h=6.626e-34; eV=1.783e-36; Msun=1.989e30
nu=lambda y: 1/(-np.expm1(-np.sqrt(y)))
meas=lambda gb: 10**(0.52*np.log10(gb)-4.19)
wb=0.0225
print(f"a0 = cH0/6 = {a0:.3e}")
print("g_bar     needed fermion mass / baryon mass   -> needed ω_f")
fs=[]
for gb in [1e-12,3e-12,1e-11,3e-11,1e-10,3e-10,1e-9]:
    f=brentq(lambda f: gb*(1+f)*nu(gb*(1+f)/a0)-meas(gb), -0.99, 1e4); fs.append(f)
    print(f"{gb:.0e}    {f:6.2f}                               {f*wb:.3f}")
# capacity: using NFW cluster densities/speeds from the phase test (inline)
rhoc=3*(70e3/Mpc)**2/(8*np.pi*G)
def nfw(M200,cc,r):
    r200=(3*M200*Msun/(4*np.pi*200*rhoc))**(1/3); rs=r200/cc; dc=200/3*cc**3/(np.log(1+cc)-cc/(1+cc)); x=r/rs
    return rhoc*dc/(x*(1+x)**2), np.sqrt(G*4*np.pi*rhoc*dc*rs**3*(np.log(1+x)-x/(1+x))/r)
def rhomax(m_eV, sig): m=m_eV*eV; return 2*m**4*(2*np.pi*sig**2)**1.5/h**3
print("\nPauli capacity vs the dark density a cluster/galaxy needs (m = 2 eV):")
for lab,(rho,v) in [("cluster 1e15, r=0.2 Mpc",nfw(1.2e15,3.7,0.2*Mpc)),("cluster 1e15, r=0.5 Mpc",nfw(1.2e15,3.7,0.5*Mpc)),("cluster 1e15, r=1.0 Mpc",nfw(1.2e15,3.7,1.0*Mpc)),
                    ("group 1e13, r=0.1 Mpc",nfw(1e13,7,0.1*Mpc)),("galaxy halo (KiDS) r=0.3 Mpc, v 220",(9e-25,220e3))]:
    cap=rhomax(2.0, v/np.sqrt(2))
    print(f"   {lab:36s} capacity/needed = {cap/rho:.2f}")
print("   SPARC galaxy outskirts (~11 kpc): m_crit ≥ 6 eV -> 2 eV capacity ≤ (2/6)^4 ≈ 1% of the galaxy's missing mass")
np.save("fermion_need.npy", np.array(fs))
