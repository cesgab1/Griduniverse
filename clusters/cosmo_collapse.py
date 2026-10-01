"""
Cosmological spherical collapse with the DBI dark fluid (weak field, physical coordinates, matter + Lambda background).
Components: baryons = collisionless shells (sorted each step, so they can cross and phase-mix); dark fluid = Lagrangian shells with
DBI pressure (artificial viscosity); fluid shells whose density runs away (> 1e4 x local mean) convert to collisionless dust shells
(the phase change: past the DBI limit the fluid is pressureless). Gravity: Newtonian + Lambda while a shell expands, MOND boost
(McGaugh nu) once it falls in (the rule used in cluster_mix.py). Start z = 50, adiabatic growing-mode top-hat whose CDM-only
collapse redshift is z_c; region simulated out to 3 Lagrangian radii.
Diagnostic at z = 0: dark mass (fluid + dust) inside the radius holding 90% of the object's baryons, relative to the cosmic share
(Omega_K/Omega_b x baryons). Want ~0 for galaxies (MOND-only), ~1 for clusters.
"""
import numpy as np, warnings, sys; warnings.filterwarnings("ignore")
G, Msun, kpc, Mpc, c, a0 = 6.674e-11, 1.989e30, 3.0857e19, 3.0857e22, 2.998e8, 1.2e-10
h = 0.68; H0 = 100*h*1e3/Mpc; Ob, OK = 0.048, 0.26; Om = Ob + OK; OL = 1 - Om; Gyr = 3.156e16
rhoc = 3*H0**2/(8*np.pi*G)
nu = lambda y: 1/(-np.expm1(-np.sqrt(np.maximum(y, 1e-30))))
def Hz(a): return H0*np.sqrt(Om/a**3 + OL)
def growth(a):   # linear growth (LCDM-like), normalised D(1)=1
    aa = np.linspace(1e-4, a, 4000); I = np.trapezoid(1/(aa*Hz(aa)/H0)**3, aa)
    a1 = np.linspace(1e-4, 1, 4000); I1 = np.trapezoid(1/(a1*Hz(a1)/H0)**3, a1)
    return Hz(a)/H0*I/(Hz(1)/H0*I1)
from scipy.integrate import quad
def t_of_a(a): return quad(lambda x: 1/(x*Hz(x)), 1e-8, a)[0]
def run(Mtot, zc, mu_inv, lam, Nb=150, Nf=250, Rfac=3.0):
    mu = 1/(mu_inv*Mpc)
    def P_of(rho):
        s = 8*np.pi*G*rho/(2*mu**2*c**2); return c**4/(8*np.pi*G)*(2*mu**2/lam)*(1 - 1/np.sqrt(1 + lam*s*s))
    ai = 1/51.; ti = t_of_a(ai); Hi = Hz(ai)
    delta_i = 1.686*growth(ai)/growth(1/(1+zc))
    rho_m_i = Om*rhoc/ai**3
    RL = (3*Mtot*Msun/(4*np.pi*rho_m_i*(1+delta_i)))**(1/3)        # physical radius of the top-hat at z_i
    Rmax = Rfac*RL
    def init(N, Omx):
        edges = np.linspace(0, Rmax, N+1); mid = edges[1:]
        dens = np.where(mid <= RL, 1+delta_i, 1.0)*Omx*rhoc/ai**3
        m = 4*np.pi/3*dens*np.diff(edges**3)
        Menc_tot = np.cumsum(4*np.pi/3*np.where(mid <= RL, 1+delta_i, 1.0)*Om*rhoc/ai**3*np.diff(edges**3))
        dl = Menc_tot/(4*np.pi/3*Om*rhoc/ai**3*mid**3) - 1
        v = Hi*mid*(1 - dl/3)
        return mid.copy(), v, m
    rb, vb, mb = init(Nb, Ob); rf, vf, mf = init(Nf, OK)
    rd = np.zeros(0); vd = np.zeros(0); md = np.zeros(0)       # dust (converted fluid)
    soft = 2*kpc*(Mtot/1e12)**(1/3)
    t = ti; tend = t_of_a(1.0)
    while t < tend:
        a_now = None
        # enclosed mass at any radius: all components
        R_all = np.concatenate([rb, rf, rd]); M_all = np.concatenate([mb, mf, md]); o = np.argsort(R_all)
        Rs, Mc = R_all[o], np.cumsum(M_all[o])
        def Menc(r):
            i = np.searchsorted(Rs, r, side="left"); return np.where(i > 0, Mc[np.maximum(i-1, 0)], 0.0)
        def grav(r, v, mself):
            M = Menc(r) + 0.5*mself
            gN = G*M/(r**2 + soft**2)
            g = np.where(v < 0, nu(gN/a0)*gN, gN)
            return -g + H0**2*OL*r
        ab = grav(rb, vb, mb)
        ad = grav(rd, vd, md) if len(rd) else np.zeros(0)
        # fluid: gravity + pressure (Lagrangian, ordered)
        af = grav(rf, vf, mf)
        edges = np.concatenate([[rf[0]*0 if len(rd) == 0 else 0.0], rf]); vol = 4*np.pi/3*np.diff(np.concatenate([[0], rf])**3)
        bad = vol <= 0
        rho = mf/np.where(bad, 1e-300, vol)
        P = P_of(rho); dv = np.diff(np.concatenate([[0], vf])); q = np.where(dv < 0, 2.0*rho*dv**2, 0.0)
        Pt = P + q; Pout = np.concatenate([Pt[1:], [Pt[-1]]]); medge = 0.5*(mf + np.concatenate([mf[1:], [mf[-1]]]))
        af = af - 4*np.pi*rf**2*(Pout - Pt)/medge
        # time step
        a_now = np.interp(t, [ti, tend], [ai, 1.0])
        drf = np.diff(np.concatenate([[0], rf])); cs = np.sqrt(np.maximum(np.diff(np.concatenate([[0], P]))/np.maximum(np.diff(np.concatenate([[0], rho])), 1e-300), 0))
        dt = 0.25*np.min(np.abs(drf)/(cs + np.abs(vf) + 1e3))
        tdyn = 0.02*np.min(np.sqrt((np.concatenate([rb, rd]) + soft)**3/(G*Mtot*Msun)))
        dt = min(dt, tdyn, 0.01*Gyr)
        vb += ab*dt; rb = np.abs(rb + vb*dt); vf += af*dt; rf = rf + vf*dt
        if len(rd): vd += ad*dt; rd = np.abs(rd + vd*dt)
        t += dt
        # conversion: fluid shells with runaway density or crossing -> dust (from the centre outward)
        vol = 4*np.pi/3*np.diff(np.concatenate([[0], rf])**3); rho = mf/np.where(vol > 0, vol, 1e-300)
        rho_mean = OK*rhoc*(Hz(1)/H0)**0 / (np.interp(t, [ti, tend], [ai, 1.0])**3)
        conv = (vol <= 0) | (rho > 1e4*rho_mean) | (rf <= 0)
        if conv.any():
            k = np.max(np.where(conv)[0]) + 1                        # convert innermost k shells
            rd = np.concatenate([rd, np.abs(rf[:k])]); vd = np.concatenate([vd, vf[:k]]); md = np.concatenate([md, mf[:k]])
            rf, vf, mf = rf[k:], vf[k:], mf[k:]
            if len(rf) < 5: break
    Mb_obj = np.sum(mb[np.argsort(rb)][:int(Nb/Rfac**3)+1]); rbs = np.sort(rb); cumb = np.cumsum(mb[np.argsort(rb)])
    r90 = np.interp(0.9*Mb_obj, cumb, rbs)
    Mdark = np.sum(mf[rf < r90]) + np.sum(md[rd < r90]); share = 0.9*Mb_obj*OK/Ob
    return r90/kpc, Mdark/share, np.sum(md)/(np.sum(md)+np.sum(mf)), (np.sum(md[rd < r90]))/share
cases = [(1e12, 2.0, "galaxy"), (1e13, 1.0, "group"), (1e14, 0.6, "cluster"), (1e15, 0.3, "massive cluster")]
mu_inv, lam = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0, float(sys.argv[2]) if len(sys.argv) > 2 else 9.2e6
print(f"1/mu = {mu_inv} Mpc, lambda_D = {lam:.1e}")
for M, zc, lab in cases:
    r90, frac, fconv, fdust = run(M, zc, mu_inv, lam)
    print(f"   {lab:16s} M = {M:.0e} Msun (z_c {zc}): baryon r90 = {r90:7.0f} kpc | dark mass inside / cosmic share = {frac:5.2f} "
          f"(of which converted dust {fdust:5.2f}) | fraction of all simulated fluid converted {fconv:.2f}", flush=True)
