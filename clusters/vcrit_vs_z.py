# Critical well depth for the DBI phase change as a function of redshift: v_crit^2 = c^2 (x_lim - x_bg(a)),
# x_bg = s/sqrt(1+lam s^2), s = s0/a^3 (background fluid gets denser in the past -> closer to the limit)
import numpy as np
G, c, Mpc = 6.674e-11, 2.998e8, 3.0857e22; H0 = 68e3/Mpc; rho = 0.26*3*H0**2/(8*np.pi*G)
for mu_inv, lam in ((300.0, 9.2e6), (22.3, 1e10)):
    mu = 1/(mu_inv*Mpc); s0 = 8*np.pi*G*rho/(2*mu**2*c**2); xl = 1/np.sqrt(lam)
    row = []
    for z in (0, 0.5, 1, 2, 3, 6):
        s = s0*(1+z)**3; x = s/np.sqrt(1+lam*s*s); row.append(f"z={z}: {c*np.sqrt(max(xl-x, 0))/1e3:6.0f}")
    print(f"1/mu={mu_inv} Mpc, lam={lam:.1e}: v_crit (km/s)  " + "  ".join(row))
print("Typical potential depths: Milky-Way-like ~500 km/s, groups ~1000, clusters ~2000+ (our MOND wells).")
