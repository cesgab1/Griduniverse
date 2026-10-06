"""Iteration 112: Coma broken down to known quantum constituents; all energy counted as mass (PREREG_112.md)."""
import os, numpy as np
c, G, h, hbar, kB, eV = 2.998e8, 6.674e-11, 6.626e-34, 1.0546e-34, 1.381e-23, 1.602e-19
MS, MPC, GYR = 1.989e30, 3.086e22, 3.156e16
mp, mn, me, mu = 1.6726e-27, 1.6749e-27, 9.109e-31, 1.6605e-27
M_need, M_gas, M_st = 1.6e15, 2.5e14, 0.5e14
X, Y = 0.75, 0.25
def run(Rmpc):
    R = Rmpc*MPC; V = 4/3*np.pi*R**3; rows = []
    Mvis = (M_gas + M_st)*MS
    # 1 particle counting: hydrogen atoms (p + e), helium-4 (2p 2n 2e)
    NH = X*Mvis/(mp + me); NHe = Y*Mvis/(4.002602*mu)
    Np, Nn, Ne = NH + 2*NHe, 2*NHe, NH + 2*NHe
    raw = Np*mp + Nn*mn + Ne*me
    bindHe = NHe*28.296e6*eV/c**2
    rows += [("1a protons + neutrons + electrons, counted one by one", raw),
             ("1b  helium nuclear binding (mass lost when He formed)", -bindHe),
             ("1c  of which quark rest masses (~1% of nucleons)", 0.01*(Np*mp + Nn*mn), "info")]
    # 2 heat of gas, mu = 0.6
    Ngas = M_gas*MS/(0.6*mp); rows.append(("2 heat of the hot gas (8.2 keV)", 1.5*Ngas*8.2e3*eV/c**2))
    # 3 gravitational binding of visible matter (uniform sphere, generous: whole visible mass)
    rows.append(("3 gravitational binding of visible matter", -0.6*G*Mvis**2/R/c**2))
    # 4 starlight: M/L ~ 5 -> L = 1e13 Lsun
    L = M_st/5*3.828e26
    rows.append(("4a starlight inside the cluster now (L x R/c)", L*R/c/c**2))
    rows.append(("4b ALL starlight ever made there (L x 12 Gyr, upper bound)", L*12*GYR/c**2))
    # 5 CMB photons in volume: u = a T^4
    rows.append(("5 relic light (CMB) inside the volume", 7.566e-16*2.7255**4*V/c**2))
    # 6 magnetic field, 5 microgauss everywhere (generous)
    rows.append(("6 magnetic field (5 microgauss throughout)", (5e-10)**2/(2*4e-7*np.pi)*V/c**2))
    # 7 cosmic rays <= 10% of heat
    rows.append(("7 cosmic rays (<= 10% of the heat)", 0.1*rows[3][1]))
    # 8 neutrinos: per mass eigenstate n_mean = 112 /cm^3 (nu + nubar); trapped max from phase space
    def nu(masses):
        tot = 0
        for m in masses:
            if m == 0: continue
            mk = m*eV/c**2; sig = 1.0e6
            n_trap = 2*0.5*(2*np.pi)**1.5*(mk*sig)**3/h**3          # g = 2, peak occupancy 1/2 (Fermi-Dirac, conserved)
            tot += (112e6 + n_trap)*V*(mk + 3.15*8.617e-5*1.945*eV/c**2*0)   # rest mass dominates
        return tot
    rows.append(("8a relic neutrinos, masses 0 / 0.0086 / 0.050 eV", nu([0, 0.0086, 0.050])))
    rows.append(("8b relic neutrinos, extreme 3 x 0.1 eV (excluded by DESI)", nu([0.1]*3), "info"))
    # 9 black holes
    rows.append(("9 central black holes (0.2% of stars)", 0.002*M_st*MS))
    # 10 measured vacuum energy, rho_L = 0.69 rho_crit (h 0.68); gravitating combination rho + 3p = -2 rho
    rhoL = 0.69*3*(68e3/MPC)**2/(8*np.pi*G)
    rows.append(("10 measured vacuum energy (dark energy): pulls as -2 x rho V", -2*rhoL*V))
    # 11 naive vacuum energy with Planck cutoff
    rhoP = c**5/(hbar*G**2)
    naive = rhoP*V
    tot = sum(r[1] for r in rows if len(r) == 2)
    out = [f"Coma inside R = {Rmpc} Mpc -- needed mass {M_need:.2e} Msun; particles: {Np:.2e} protons, {Nn:.2e} neutrons, {Ne:.2e} electrons", ""]
    for r in rows:
        tag = " (not added)" if len(r) == 3 else ""
        out.append(f"  {r[0]:62s} {r[1]/MS:+.3e} Msun  = {r[1]/MS/M_need:+.2e} of needed{tag}")
    beyond = tot/MS - raw/MS
    out += ["", f"  TOTAL of all known constituents: {tot/MS:.3e} Msun = {tot/MS/M_need:.3f} of the needed mass  -> short x{M_need/(tot/MS):.2f}",
            f"  everything beyond plain protons/neutrons/electrons: {beyond:+.2e} Msun = {beyond/M_need:+.1e} of needed",
            f"  11 naive quantum vacuum energy (Planck cutoff) in the volume: {naive/MS:.1e} Msun = {naive/MS/M_need:.0e} x needed (and repulsive)"]
    return out
here = os.path.dirname(os.path.abspath(__file__)); out = ["Iteration 112 -- Coma, broken down to its known quantum constituents"]
for R in (2.5, 2.0, 3.0): out += [""] + run(R) if R == 2.5 else ["", f"[R = {R} Mpc] " + run(R)[-3].strip()]
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter112_inventory.txt"), "w").write(txt + "\n")
