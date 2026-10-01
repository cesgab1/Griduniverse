"""
Does Khronon's own dynamics make the fluid around a galaxy settle into exactly the MOND halo?
Equations (Blanchet & Skordis 2024, arXiv:2404.06584, weak field; checked independently, see README):
  Delta phi = 4 pi G (rho_m + rho_tau)                                   (3.11)
  rho_tau = -(1/4 pi G) [ div(J_Y grad Xi) + mu^2 Xi ]                   (3.13)  dust + MOND 'phantom' in ONE density
  d rho_tau/dt + div(rho_tau v) = 0                                      (3.16)  conserved, carried by the khronon flow
  Dv/Dt = -grad phi + grad Xi                                            (from Xi = phi - sigma_dot + |grad sigma|^2/2, v = -grad sigma)
  f = 1 + J_Y is the MOND function of |grad Xi|/a0 (we use the inverse of nu(y) = 1/(1 - exp(-sqrt y)), as in all our fits).
Spherical, galaxy scales (mu^2 term negligible for 1/mu >~ 20 Mpc; it carries the uniform background):
  (1 - f(|Xi'|/a0)) Xi' = G dM_tau(<r)/r^2,  dM_tau = fluid mass in excess of the cosmic mean.
Static equilibrium (v = 0) needs Xi' = phi'  =>  f(g/a0) g = G M_b/r^2: exactly MOND, with the fluid excess = MOND phantom mass.
Run: Lagrangian fluid shells in an expanding universe (start z = 30, Hubble flow) around a galaxy of baryons M_b (fixed, plus
the mean baryons of each shell's region beyond), to z = 0. Compare: (a) Khronon dynamics; (b) the model used in our earlier
tests (fluid as dust with its own Newtonian gravity, baryons pulling with MOND) -- the double-counting suspect.
Units kpc, Msun, Gyr.
"""
import numpy as np, sys
from scipy.integrate import quad
k = 1.0227; G = 4.30e-6*k**2; a0 = 1.2e-10*3.086e19/1e6*k**2
h = 0.68; H0 = 0.1*h*k; Om, Ob, OL = 0.31, 0.049, 0.69; Oc = Om - Ob
rhoc = 3*H0**2/(8*np.pi*G)
nu = lambda y: 1/(1 - np.exp(-np.sqrt(np.maximum(y, 1e-14))))
ys = np.logspace(-8, 4, 4000); xs = ys*nu(ys); fs = 1/nu(ys)              # g/a0 = x, f(x) = g_N/g
f_of = lambda x: np.interp(x, xs, fs)
hx = (1 - fs)*xs; ipk = np.argmax(hx); hmax = hx[ipk]                      # (1-f) g peaks at ~0.65 a0
def Xi_prime(s):                                                         # solve (1-f(g/a0)) g = s on the small-g branch
    sg = np.sign(s); s = np.abs(s)/a0
    x = np.interp(np.minimum(s, hmax), hx[:ipk+1], xs[:ipk+1])
    return sg*x*a0, s > hmax
def tz(z): return quad(lambda a: 1/(a*H0*np.sqrt(Om/a**3 + OL)), 0, 1/(1+z))[0]
def run(Mb, mode, N=300, zi=30.0, Rmax=4000.0):
    ai = 1/(1+zi); R = np.geomspace(5, Rmax, N)                          # comoving radii (kpc)
    edges = np.concatenate([[0], R]); dV = 4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
    m_f = Oc*rhoc*dV                                                     # fluid mass per shell
    Mb_enc = lambda Rc: Mb + Ob*rhoc*4/3*np.pi*np.maximum(Rc**3 - (3*Mb/(4*np.pi*Ob*rhoc))**1, 0)
    Mb_shell = np.maximum(Mb_enc(R), Mb)                                 # baryons inside each shell's comoving radius (follow shells)
    r = ai*R; v = H0*np.sqrt(Om/ai**3 + OL)*r
    t = tz(zi); out = {}; eps = 3.0; broke = 0
    marks = {tz(0.9): "z=0.9", tz(0.0): "z=0"}
    while True:
        a_now = None
        order = np.argsort(r); Mf = np.empty(N); Mf[order] = np.cumsum(m_f[order])
        Mb_r = np.interp(r, ai*R*0 + r[order], Mb_shell[order]) if False else Mb_shell     # baryons ride with shells
        # mean cosmic fluid inside radius r at time t: need a(t); integrate a with Friedmann alongside
        rr = np.sqrt(r**2 + eps**2)
        gN_b = G*Mb_r/rr**2
        if mode == "khronon":
            dM = Mf - 4/3*np.pi*r**3*Oc*rhoc/aa(t)**3
            Xp, bad = Xi_prime(G*dM/rr**2); broke = max(broke, bad.sum())
            acc = -G*(Mb_r + Mf)/rr**2 + Xp + OL*H0**2*r
        else:                                                            # old model: baryons MOND, fluid Newtonian dust
            acc = -(nu(gN_b/a0)*gN_b) - G*Mf/rr**2 + OL*H0**2*r
        dt = max(min(0.005, 0.05*np.min(rr/(np.abs(v) + 1e-3))), 2e-4)
        v += acc*dt; r += v*dt; t += dt
        neg = r < 2.0; r[neg] = 2.0; v[neg] = np.abs(v[neg])
        for tm, lab in list(marks.items()):
            if t >= tm and lab not in out:
                order = np.argsort(r); Mfs = np.cumsum(m_f[order]); rs = r[order]
                out[lab] = (rs.copy(), Mfs.copy(), aa(t))
        if t >= tz(0.0): break
    return out, broke
_ag = np.linspace(1e-3, 1.0, 20000); _tg = np.array([tz(1/a - 1) for a in _ag[::200]])
def aa(t): return np.interp(t, _tg, _ag[::200])
if __name__ == "__main__":
    for Mb in (1e10, 1e11):
        res = {}
        for mode in ("khronon", "old"):
            res[mode], br = run(Mb, mode)
            if mode == "khronon": print(f"\nMb = {Mb:.0e} Msun   (khronon: shells beyond the 0.65 a0 cap at some step: {br})")
        print("   r (kpc) | MOND phantom | Khronon excess fluid z=0.9 / z=0 | old model excess z=0.9 / z=0")
        for rk in (10, 30, 100, 300, 1000):
            gN = G*Mb/rk**2; ph = (nu(gN/a0) - 1)*Mb
            row = []
            for mode in ("khronon", "old"):
                for lab in ("z=0.9", "z=0"):
                    rs, Mfs, av = res[mode][lab]
                    exc = np.interp(rk, rs, Mfs) - 4/3*np.pi*rk**3*Oc*rhoc/av**3
                    row.append(exc)
            print(f"   {rk:7d} | {ph:12.2e} | {row[0]:11.2e} / {row[1]:9.2e}     | {row[2]:9.2e} / {row[3]:9.2e}")
        sys.stdout.flush()
