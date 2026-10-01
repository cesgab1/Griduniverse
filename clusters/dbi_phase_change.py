"""
DBI Khronon fluid in static galaxies vs clusters (weak field). Bernoulli for the khronon flow: Q - 1 = x = x0 + |phi|/c^2 (phi -> 0 far away).
DBI: K' = 2 mu^2 x / sqrt(1 - lam x^2) diverges at x_lim = 1/sqrt(lam): fluid density -> infinity, pressure saturates -> pressureless dust
that can collapse ("phase change"). Critical well depth: |phi_c| = c^2 (x_lim - x0) -> v_crit = sqrt(|phi_c|).
Background today: x0 = s0/sqrt(1 + lam s0^2), s0 = I0/(2 mu^2), I0 = 8 pi G rho_K (CLASS units 3 rho).
Potentials: MOND (McGaugh nu) from baryons, zero at R = 3 Mpc, for a Milky-Way-like galaxy, a group and two clusters.
"""
import numpy as np, warnings; warnings.filterwarnings("ignore")
exec(open("khronon_mass_term.py").read().split('print("Cluster')[0])
c = 2.998e8; I0 = 3*0.12/0.68**2*(68/2.998e5)**2; mu = 1/22.3
def scaled(Mgas, rc_kpc, Mstar, a_kpc):
    def Mb(r):
        rc = rc_kpc*kpc; x = r/rc; m = x - np.arctan(x); m1 = 1e3*kpc/rc - np.arctan(1e3*kpc/rc)
        return Mgas*Msun*m/m1 + Mstar*Msun*r**2/(r + a_kpc*kpc)**2
    return Mb
objs = [("Milky-Way-like galaxy", scaled(1e10, 100, 6e10, 3)), ("massive spiral", scaled(3e10, 150, 2e11, 5)),
        ("galaxy group", scaled(5e12, 100, 1e12, 10)), ("cluster 1e14 Msun gas", scaled(1e14, 150, 1e12, 20)),
        ("massive cluster 3e14 gas", scaled(3e14, 250, 3e12, 30))]
r = np.geomspace(1*kpc, 3*Mpc, 3000)
pots = {}
for lab, Mb in objs:
    g, gN, _ = solve(Mb, None, 3*Mpc, r)
    phi = -np.concatenate([np.cumsum((g[1:]*np.diff(r))[::-1])[::-1], [0]])
    vflat = np.sqrt(np.interp(50*kpc if "galaxy" in lab or "spiral" in lab else 500*kpc, r, g*r))/1e3
    pots[lab] = (np.abs(phi[0]), vflat)
print(f"{'lambda_D':>9s} {'v_crit':>8s} | " + " | ".join(f"{l[:22]:>22s}" for l, _ in objs))
print(f"{'':>9s} {'(km/s)':>8s} | " + " | ".join(f"{'sqrt|phi_0| ' + str(int(np.sqrt(pots[l][0])/1e3)) + ' km/s':>22s}" for l, _ in objs))
for lam in (3e9, 1e10, 3e10):
    s0 = I0/(2*mu**2); x0 = s0/np.sqrt(1+lam*s0**2); xl = 1/np.sqrt(lam); vc = c*np.sqrt(xl - x0)/1e3
    def rho_ratio(phi_abs):
        x = x0 + phi_abs/c**2
        if x >= xl: return np.inf
        Kp = lambda xx: 2*mu**2*xx/np.sqrt(1 - lam*xx**2)
        return Kp(x)*(1+x)/(Kp(x0)*(1+x0))       # rho ~ Q K' (K negligible)
    cells = []
    for l, _ in objs:
        rr = rho_ratio(pots[l][0])
        cells.append("COLLAPSES (dust)" if not np.isfinite(rr) else f"fluid x{rr:.2f} (no halo)")
    print(f"{lam:9.0e} {vc:8.0f} | " + " | ".join(f"{cc:>22s}" for cc in cells))
