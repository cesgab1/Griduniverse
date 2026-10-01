"""
Planck stars (Rovelli & Vidotto 2014; Haggard & Rovelli 2015) in grid terms: collapse stops when the grid reaches Planck DENSITY
(links cannot be squeezed further) at radius r ~ (M/m_P)^(1/3) l_P, bounces, and the outside sees a long-lived black hole.
Lifetimes: bounce tau = 4 k (M/m_P)^2 t_P (k ~ 0.05, Christodoulou et al 2014; range up to ~1e22 scanned by Barrau et al.)
           Hawking  tau = 5120 pi G^2 M^3 / (hbar c^4)   (the M^3 scenario: evaporate first, then a Planck-mass white-hole remnant)
"""
import numpy as np
G, hbar, c = 6.674e-11, 1.0546e-34, 2.998e8
mP, tP, lP = np.sqrt(hbar*c/G), np.sqrt(hbar*G/c**5), np.sqrt(hbar*G/c**3)
Msun, yr, age = 1.989e30, 3.156e7, 13.8e9
objs = [("primordial hole exploding today (k=0.05)", None), ("stellar black hole, 10 Msun", 10*Msun),
        ("Sgr A* (Milky Way centre), 4.3e6 Msun", 4.3e6*Msun), ("M87*, 6.5e9 Msun", 6.5e9*Msun), ("largest known, ~4e10 Msun", 4e10*Msun)]
Mtoday = mP*np.sqrt(age*yr/(4*0.05*tP))
print(f"Planck length {lP:.2e} m, Planck mass {mP:.2e} kg\n")
print(f"{'object':44s} {'mass kg':>10s} {'core radius':>12s} {'bounce tau (M^2)':>17s} {'Hawking tau (M^3)':>18s}")
for lab, M in objs:
    M = Mtoday if M is None else M
    r = (M/mP)**(1/3)*lP
    tb = 4*0.05*(M/mP)**2*tP/yr; th = 5120*np.pi*G**2*M**3/(hbar*c**4)/yr
    print(f"{lab:44s} {M:10.2e} {r:10.2e} m {tb:14.2e} yr {th:15.2e} yr")
print(f"\nAge of universe {age:.1e} yr. Every astrophysical black hole, including all supermassive ones, is far from bouncing in either scenario.")
