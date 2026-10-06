"""Drill D4 (exploratory): Coma's hot-gas (hydrostatic) mass vs weak-lensing masses, and what that does to the
emergent-gravity ratio of iteration 113 (its disagreement with Tamosiunas+2019)."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
G, MS, KPC, MPC, eV, mp = 6.674e-11, 1.989e30, 3.086e19, 3.086e22, 1.602e-19, 1.6726e-27
h = 0.7; rhoc = 3*(100*h*1e3/MPC)**2/(8*np.pi*G)
beta, rc, ne0, kT = 0.75, 291*KPC, 3.44e-3*1e6, 8.2e3*eV
r = np.linspace(1, 4000, 40000)*KPC
Mhse = 3*beta*kT*r**3/(G*0.6*mp*(r**2 + rc**2))
Mgas = cumulative_trapezoid(4*np.pi*r**2*1.17*mp*ne0*(1 + (r/rc)**2)**(-1.5*beta), r, initial=0)
MB = Mgas + 5e13*MS*Mgas/np.interp(2500*KPC, r, Mgas)
MD = np.sqrt(np.clip(2.998e8*(70e3/MPC)*r**2/(6*G)*np.gradient(MB*r, r), 0, None)); Meg = MB + MD
def at_delta(Mprof, D):
    i = np.argmin(np.abs(Mprof/(4/3*np.pi*r**3) - D*rhoc)); return r[i], Mprof[i]
out = ["Drill 4 -- Coma: hot-gas mass vs light-bending mass"]
for lab, D, M, lo, hi in [("Kubo+2007 SDSS, M200", 200, 1.88e15/h, 0.56e15/h, 0.65e15/h), ("Okabe+2014 Subaru, Mvir (~100x critical)", 100, 8.42e14/h, 2.42e14/h, 4.17e14/h)]:
    rl = (3*M*MS/(4*np.pi*D*rhoc))**(1/3)               # radius implied by the lensing mass itself
    mh = np.interp(rl, r, Mhse)/MS; me = np.interp(rl, r, Meg)/MS
    out.append(f"  {lab}: {M:.2e} (-{lo:.1e} +{hi:.1e}) Msun inside {rl/MPC/1e-3/1e3:.2f} Mpc;  hot-gas mass there {mh:.2e} (lensing/hot-gas x{M/mh:.2f});"
               f"  emergent gravity there {me:.2e} -> EG/lensing x{me/M:.2f}")
txt = "\n".join(out); print(txt); open("drill4_coma_lensing.txt", "w").write(txt + "\n")
