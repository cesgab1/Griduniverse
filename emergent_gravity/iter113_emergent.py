"""Iteration 113: Verlinde emergent gravity in Coma + four entropies (PREREG_113.md)."""
import os, numpy as np
from scipy.integrate import cumulative_trapezoid
c, G, h, hbar, kB, eV = 2.998e8, 6.674e-11, 6.626e-34, 1.0546e-34, 1.381e-23, 1.602e-19
MS, KPC, mp, me = 1.989e30, 3.086e19, 1.6726e-27, 9.109e-31
H0 = 70e3/(1e3*KPC); L = c/H0; lP2 = hbar*G/c**3
beta, rc, ne0, kT = 0.75, 291*KPC, 3.44e-3*1e6, 8.2e3*eV
r = np.linspace(1, 3000, 30000)*KPC
ne = ne0*(1 + (r/rc)**2)**(-1.5*beta)
Mgas = cumulative_trapezoid(4*np.pi*r**2*1.17*mp*ne, r, initial=0)
def stars(kind):
    if kind == "gas-shape": s = Mgas/np.interp(2500*KPC, r, Mgas)
    else: a = 150*KPC; s = (r/(r + a))**2/((2500*KPC)/(2500*KPC + a))**2
    return 5e13*MS*s
Mmeas = 3*beta*kT*r**3/(G*0.6*mp*(r**2 + rc**2))
radii = [100, 300, 500, 1000, 1500, 2500]
out = ["Iteration 113 -- emergent gravity (Verlinde) in Coma; no free parameter (H0 = 70)",
       f"check: gas mass within 2.5 Mpc {np.interp(2500*KPC, r, Mgas)/MS:.2e} Msun (project 2.5e14); hydrostatic mass at 2.5 Mpc "
       f"{np.interp(2500*KPC, r, Mmeas)/MS:.2e} (project 1.6e15); Verlinde a0 = cH0/6 = {c*H0/6:.2e} m/s^2", ""]
for kind in ("gas-shape", "hernquist"):
    MB = Mgas + stars(kind)
    MD = np.sqrt(np.clip(c*H0*r**2/(6*G)*np.gradient(MB*r, r), 0, None))
    out.append(f"stars {kind}:   r [kpc] | visible | emergent extra | EG total | measured (hydrostatic) | EG/measured")
    for R in radii:
        i = np.searchsorted(r, R*KPC); mb, md, mm = MB[i]/MS, MD[i]/MS, Mmeas[i]/MS
        out.append(f"   {R:5d} | {mb:.2e} | {md:.2e} | {mb+md:.2e} | {mm:.2e} | {(mb+md)/mm:.2f}")
    out.append("")
# ---- entropy ledger at R = 2.5 Mpc ----
R = 2500*KPC; i = np.searchsorted(r, R); MB = Mgas + stars("gas-shape")
def st(n, m):           # Sackur-Tetrode entropy density / k
    lam = h/np.sqrt(2*np.pi*m*kT); return n*(np.log(1/(n*lam**3)) + 2.5)
sgas = st(ne, me) + st(0.83*ne, mp) + st(0.083*ne, 4*mp)
S_gas = np.trapezoid((4*np.pi*r**2*sgas)[:i], r[:i])
N_star_baryons = 5e13*MS/mp; S_star = 15*N_star_baryons
V = 4/3*np.pi*R**3; s_cmb = 4/3*7.566e-16*2.7255**3/kB; S_cmb = s_cmb*V; S_nu = 0.9545*S_cmb
S_bh1 = 4*np.pi*G*(2e10*MS)**2/(hbar*c); S_bh2 = 1000*4*np.pi*G*(1e8*MS)**2/(hbar*c)
A = 4*np.pi*R**2; S_area = A/(4*lP2); S_DE = S_area*R/L; S_M = 2*np.pi*MB[i]*c*R/hbar
out += ["ENTROPY LEDGER for Coma inside 2.5 Mpc (units of k_B)",
        "  thermodynamic:",
        f"    hot gas (Sackur-Tetrode, electrons + ions)          {S_gas:.1e}",
        f"    stars (~15 k per baryon)                             {S_star:.1e}",
        f"    relic light (CMB) in the volume                      {S_cmb:.1e}",
        f"    relic neutrinos in the volume                        {S_nu:.1e}",
        "  coarse-grained gravitational (black holes, Bekenstein-Hawking):",
        f"    one 2e10 Msun black hole (NGC 4889 class)            {S_bh1:.1e}",
        f"    1000 x 1e8 Msun black holes                          {S_bh2:.1e}",
        "  entanglement / holographic:",
        f"    area law A/4l_P^2 (= maximum the region can hold)    {S_area:.1e}",
        f"    de Sitter volume-law part S_DE = area x r/L          {S_DE:.1e}",
        f"    entropy displaced by Coma's visible matter S_M       {S_M:.1e}",
        "  von Neumann entropy of the whole (cluster + rest of universe, pure state): 0, and constant in time",
        ""]
gN = G*MB/r**2
out.append("Where emergent gravity switches on: S_M/S_DE = g_N/(cH0/2)  (< 1 -> 'weak pull', apparent dark matter appears)")
out.append("   " + ", ".join(f"{R} kpc: {np.interp(R*KPC, r, gN)/(c*H0/2):.3f}" for R in radii))
txt = "\n".join(out); print(txt); open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "iter113_emergent.txt"), "w").write(txt + "\n")
