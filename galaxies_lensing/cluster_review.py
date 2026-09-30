# The cluster weak spot, quantified, and grid remedies checked
import numpy as np
c=2.998e8; Mpc=3.0857e22; H0=67.7e3/Mpc; G=6.674e-11; Msun=1.989e30; kpc=3.0857e19
a0=c*H0/6
nu=lambda y: 1/(1-np.exp(-np.sqrt(y)))
clash=lambda gb: 10**(0.52*np.log10(gb)-4.19)     # BCG+cluster RAR (Tian+2024), 4.9% scatter
print(f"slack rule a0 = cH0/6 = {a0:.2e};  clusters behave as if a0 = 2.0e-9 (×{2e-9/a0:.0f})\n")
print("g_bar       measured (clusters)   slack rule    missing factor in pull   ...in mass")
for gb in [1e-12,3e-12,1e-11,3e-11,1e-10,3e-10,1e-9]:
    m=clash(gb); p=gb*nu(gb/a0)
    print(f"{gb:.0e}      {m:.2e}            {p:.2e}       ×{m/p:4.1f}                   ×{m/p:4.1f}")
# Remedy 1: extra baryons. Deep regime g ∝ sqrt(g_bar): need baryons ×(factor)^2
f=clash(1e-11)/(1e-11*nu(1e-11/a0))
print(f"\n1. hidden baryons: need ×{f**2:.0f} more baryons than seen (clusters already hold the cosmic baryon share ~0.16 of total ΛCDM mass; ×{f**2:.0f} would exceed all baryons in the universe per unit mass)")
# Remedy 2: neutrinos (Sanders/Angus route) vs DESI+CMB
print("2. neutrinos: MOND clusters need ~2 eV neutrinos each (Sanders 2003) or an 11 eV sterile one (Angus 2009); DESI DR2 + CMB: sum < 0.064 eV -> active route excluded")
# Remedy 3: each galaxy keeps its own slack and they add linearly (your 'two dimples' picture) - Coma-like cluster
Ngal, Mstar, Mgas = 1000, 1e10*Msun, 1e14*Msun
Mtot=Ngal*Mstar+Mgas
ph_lin = Ngal*np.sqrt(G*Mstar*a0)/G    # phantom mass per unit radius (kg/m), linear sum over galaxies
ph_mond = np.sqrt(G*Mtot*a0)/G
print(f"3. per-galaxy slack adding linearly (Coma-like: 1000 galaxies of 1e10 Msun + 1e14 Msun gas): phantom mass ×{(ph_lin+np.sqrt(G*Mgas*a0)/G)/ph_mond:.0f} vs MOND"
      f" -> overshoots the ×{f:.1f} needed, and the neighbour test (EFE, Chae+2021) showed slack combines non-linearly")
# Remedy 4: the CDM-like component our own CMB fits already contain
Oc,Ob=0.264,0.049
print(f"4. our CMB fits already contain a cold-dark-matter-like component (ω_c = 0.12, {Oc/Ob:.1f}× baryons). Clusters' total/baryon ≈ {(Oc+Ob)/Ob:.1f}, matching that share;"
      "\n   the framework has not yet said what this component is in grid terms, nor why it would gather in clusters but not dominate disc galaxies")
