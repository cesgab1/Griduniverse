# Hunting the missing factor in a0. Each refinement is a physical ingredient of the grid, applied to the earlier rules.
import numpy as np
c=2.998e8; Mpc=3.0857e22; G=6.674e-11
H0=67.7e3/Mpc; OL=0.689
rhoDE=OL*3*H0**2/(8*np.pi*G)
meas={'SPARC (z≈0)':(1.15e-10,0.1e-10),'KiDS discs (z≈0.25)':(1.07e-10,0.12e-10)}
rows=[]
g_point=c*np.sqrt(8*np.pi*G*rhoDE)                     # field energy density = tension (earlier)
rows.append(("1. pull energy density = tension, at a point (earlier)", g_point))
rows.append(("2. + mosaic: random link directions, only the share along the pull resists (1/3 of the energy)", g_point/np.sqrt(3)))
# 3. whole-volume comparison: field energy outside r vs tension energy inside r
g_vol=c*np.sqrt(8*np.pi*G*rhoDE/3)
rows.append(("3. whole volume: field energy outside r = tension inside r", g_vol))
rows.append(("4. riptide exchange counting, loops (a0 = cH)", c*H0))
rows.append(("5. elastic grid: exchange scale cH × volume/area (1/3) × elastic energy ½ = cH/6  [Verlinde 2016]", c*H0/6))
print(f"{'rule':100s} {'a0':>9s}  " + "  ".join(f"{k}" for k in meas))
for lab,g in rows:
    print(f"{lab:100s} {g:9.2e}  " + "  ".join(f"×{g/v:5.2f} ({(g-v)/e:+5.1f}σ)" for v,e in meas.values()))
# onset reading: where the fitted slack law gives only a small extra pull
nu=lambda y: 1/(1-np.exp(-np.sqrt(y)))
print(f"\nOnset reading of rule 1: at g = {g_point:.2e} the fitted slack law gives extra pull of {100*(nu(g_point/1.15e-10)-1):.1f}%  (slack 'begins' there)")
print(f"Onset reading of rule 3: at g = {g_vol:.2e} -> extra pull {100*(nu(g_vol/1.15e-10)-1):.1f}%")
