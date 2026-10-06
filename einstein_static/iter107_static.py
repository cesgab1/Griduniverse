import numpy as np
G, c, hb, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.3807e-23
lP = np.sqrt(hb*G/c**3); rhoP = c**5/(hb*G**2); TP = np.sqrt(hb*c**5/G)/kB; EP = np.sqrt(hb*c**5/G)
out = []
# Q1
R = lP; rho = 3*c**2/(16*np.pi*G*R**2)
out.append(f"Q1 Planck-radius static universe: density {rho/rhoP:.3f} Planck densities (grid cap 0.41); Lambda = {8*np.pi*G*rho/c**2*lP**2:.2f} / l_P^2")
# Q2: radiation entropy of a closed static universe of radius R: volume 2 pi^2 R^3; s = (4/3) rho c^2 / T; rho c^2 = (pi^2/30) g (kT)^4/(hbar c)^3
g = 106.75
def T_of_rho(rho): return ((rho*c**2)*(hb*c)**3/(np.pi**2/30*g))**0.25/kB
T = T_of_rho(rho); S = (4/3)*rho*c**2*2*np.pi**2*R**3/(kB*T)
S_need = 1e89
out.append(f"Q2 its temperature {T/TP:.2f} T_P, entropy {S:.1f} (units of k_B); our visible universe ~{S_need:.0e} -> must create x{S_need/S:.0e};"
           f" e-folds of conversion ~ {np.log(S_need/S)/3:.0f}")
# Q3
out.append("Q3 our law rho_DE ~ |adot|^-1/2 -> infinite at a static phase (adot = 0): incompatible unless the law is cut off near the cap")
# Q4: area-law thermal seeds. <dE^2> = T^2 C_V, C_V = (R/l)^2 ; Phi ~ G dE/(R c^4) -> P_Phi ~ (G kT/(c^4 l))^2 = ((T/T_P)(l_P/l))^2
Pzeta = 2.1e-9; PPhi = (3/5)**2*Pzeta; lcell = 1.67*lP
Tneed = np.sqrt(PPhi)*(lcell/lP)*TP
out.append(f"Q4 area-law thermal seeds: scale-invariant; amplitude needs T = {Tneed/TP:.1e} T_P = {Tneed:.1e} K (cell 1.67 l_P)")
# static Einstein universe at that temperature: radius from balance; can it hold our entropy?
rho_T = (np.pi**2/30*g*(kB*Tneed)**4/(hb*c)**3)/c**2
R_T = np.sqrt(3*c**2/(16*np.pi*G*rho_T)); S_T = (4/3)*rho_T*c**2*2*np.pi**2*R_T**3/(kB*Tneed)
# radius needed to hold our entropy at that temperature
s_dens = (4/3)*rho_T*c**2/(kB*Tneed); R_need = (S_need/(2*np.pi**2*s_dens))**(1/3)
out.append(f"   Einstein-static at that T: radius {R_T/lP:.1e} l_P, entropy {S_T:.1e}; to hold 1e89 it needs radius {R_need/lP:.1e} l_P"
           f" -> mismatch x{R_need/R_T:.0e} in size (x{S_need/S_T:.0e} in entropy)")
txt = "\n".join(out); print(txt); open("iter107_static.txt", "w").write(txt + "\n")
