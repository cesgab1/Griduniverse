"""
Time-dependent test of the DBI phase change (weak field, spherical, Lagrangian shells).
Dark fluid: irrotational barotropic fluid with Khronon DBI equation of state (enthalpy h = c^2 x, x = Q-1):
   K'(x) = 2 mu^2 x / sqrt(1 - lam x^2) = 8 pi G rho / c^2   ->  x(rho) = s/sqrt(1+lam s^2), s = 8 pi G rho/(2 mu^2 c^2)
   P = c^4 K/(8 pi G), K = (2 mu^2/lam)(1 - 1/sqrt(1+lam s^2))      (pressure saturates as rho -> infinity)
Baryons: static well (already formed), MOND (McGaugh nu) gravity. Fluid feels the total MOND field of baryons + its own excess mass.
Start: fluid at the cosmic mean density, at rest, out to 3 Mpc (outer shells held fixed = reservoir). Evolve 10 Gyr.
Outcome: 'smooth' = settles to hydrostatic profile with modest overdensity; 'collapse' = runaway infall / shell crossing.
Parameters: best point of the cosmology scan (1/mu = 300 Mpc, lambda_D = 9.2e6) plus the 22.3 Mpc / 1e10 point.
"""
import numpy as np, warnings, sys; warnings.filterwarnings("ignore")
exec(open("dbi_phase_change.py").read().split("r = np.geomspace")[0])     # objs, scaled(), solve(), nu, constants
c = 2.998e8; Gyr = 3.156e16
rho_bar = 0.26*3*(68e3/Mpc)**2/(8*np.pi*G)
def run(Mb, mu_inv, lam, N=300, Rout=3.0, tend=10.0):
    mu = 1/(mu_inv*Mpc)
    def x_of(rho):
        s = 8*np.pi*G*rho/(2*mu**2*c**2); return s/np.sqrt(1+lam*s*s), s
    def P_of(rho):
        x, s = x_of(rho); return c**4/(8*np.pi*G)*(2*mu**2/lam)*(1 - 1/np.sqrt(1+lam*s*s))
    # Lagrangian shells: edges re[0..N], mass per shell from uniform background
    re = np.linspace(0, Rout*Mpc, N+1); m = 4*np.pi/3*rho_bar*np.diff(re**3)
    r = re[1:].copy(); v = np.zeros(N); M_enc_bg = lambda rr: 4*np.pi/3*rho_bar*rr**3
    def accel(r, v):
        edges = np.concatenate([[0], r]); vol = 4*np.pi/3*np.diff(edges**3)
        if np.any(vol <= 0): return None, None
        rho = m/vol; P = P_of(rho)
        Mf = np.cumsum(m)                                   # fluid mass inside each shell edge
        gN = G*(Mb(r) + (Mf - M_enc_bg(r)))/r**2           # baryons + fluid excess over background
        g = np.sign(gN)*nu(np.abs(gN)/a0)*np.abs(gN)
        # pressure force on shell edge i: -(4 pi r^2)(P_out - P_in)/m_edge ; artificial viscosity q
        dv = np.diff(np.concatenate([[0], v])); q = np.where(dv < 0, 2.0*rho*dv**2, 0.0)
        Pt = P + q; Pout = np.concatenate([Pt[1:], [Pt[-1]]]); medge = 0.5*(m + np.concatenate([m[1:], [m[-1]]]))
        a = -g - 4*np.pi*r**2*(Pout - Pt)/medge
        return a, rho
    t, dt = 0.0, 0.0; hist = []
    rho0c = None
    while t < tend*Gyr:
        a, rho = accel(r, v)
        if a is None: return "collapse (shell crossing)", t/Gyr, hist
        x, _ = x_of(rho); cs = np.sqrt(np.maximum(np.gradient(P_of(rho))/np.maximum(np.gradient(rho), 1e-300), 0)) + 1e3
        dr = np.diff(np.concatenate([[0], r])); dt = 0.2*np.min(dr/(cs + np.abs(v) + 1e3))
        dt = min(dt, 0.02*Gyr)
        v += a*dt; v[-20:] = 0; r += v*dt; t += dt
        if len(hist) == 0 or t/Gyr - hist[-1][0] > 0.5: hist.append((t/Gyr, rho[0]/rho_bar, np.mean(rho[:10])/rho_bar))
        if rho[0]/rho_bar > 1e7: return "collapse (runaway density)", t/Gyr, hist
    return "settled", t/Gyr, hist
cases = [(300.0, 9.2e6), (22.3, 1.0e10)]
for mu_inv, lam in cases:
    I0 = 8*np.pi*G*rho_bar/c**2; s0 = I0/(2*(1/(mu_inv*Mpc))**2); x0 = s0/np.sqrt(1+lam*s0**2); vcrit = c*np.sqrt(1/np.sqrt(lam) - x0)/1e3
    print(f"\n1/mu = {mu_inv} Mpc, lambda_D = {lam:.1e}: static critical depth v_crit = {vcrit:.0f} km/s")
    for lab, Mb in objs:
        res, tt, hist = run(Mb, mu_inv, lam)
        last = hist[-1] if hist else (0, np.nan, np.nan)
        print(f"   {lab:26s}: {res:28s} at t = {tt:5.2f} Gyr | central fluid density {last[1]:.3g} x mean (inner 10 shells {last[2]:.3g})", flush=True)
