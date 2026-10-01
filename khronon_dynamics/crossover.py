"""
Crossover window with one-fluid bookkeeping (Khronon weak field, arXiv:2404.06584 eqs 3.11, 3.13, 3.16):
   div((1-f) grad Xi) - mu^2 Xi = 4 pi G drho_tau,   Dv/Dt = -grad phi + grad Xi,   f = MOND function of |grad Xi|/a0
(J_Y = f - 1). The fluid feels its own gravity on scales > 1/mu (CDM-like) and not on scales < 1/mu (MOND regime), so
crossover length lambda = 2 pi/mu. A time-dependent crossover is the only way to meet:
  Lyman-alpha (z ~ 3)  -> CDM-like above ~0.3 Mpc comoving: lambda(z=3) small
  MOND galaxies, KiDS  -> J-term regime out to ~0.3-1 Mpc around galaxies by z ~ 0.9 (lensing same at z 0.9 and 0.3)
  supply               -> fluid gathered while CDM-like must then relax to the MOND equilibrium (excess pushed out)
Schedule (comoving): lam3 at z = 3 (scaling as a^1/2 before, constant sound speed), then a power of a up to lam0 at z = 0.
Part A, linear cosmology: fluid self-gravity factor mu^2/(k^2 + mu^2) (physical k); baryons feel full gravity.
  -> sigma8 ratio and Lyman-alpha power ratio at z = 3, k = 5/Mpc, vs CDM.
Part B, one galaxy (Mb = 1e11 Msun, fixed from z = 30): Lagrangian fluid shells, Xi from the screened, nonlinear equation
  (finite volume, Picard on f), compare fluid excess with the MOND halo at 30/100/300 kpc at z = 2, 0.9, 0.
Units kpc, Msun, Gyr.
"""
import numpy as np, sys, os
DTMIN = float(os.environ.get("DTMIN", "5e-4")); NOXI = bool(os.environ.get("NOXI"))
from scipy.integrate import quad, solve_ivp
from scipy.linalg import solve_banded
k = 1.0227; G = 4.30e-6*k**2; a0 = 1.2e-10*3.086e19/1e6*k**2
h = 0.68; H0 = 0.1*h*k; Om, Ob, OL = 0.31, 0.049, 0.69; Oc = Om - Ob; rhoc = 3*H0**2/(8*np.pi*G)
nu = lambda y: 1/(1 - np.exp(-np.sqrt(np.maximum(y, 1e-14))))
ys = np.logspace(-10, 6, 5000); xs = ys*nu(ys); fs = 1/nu(ys)
f_of = lambda x: np.interp(x, xs, fs, left=fs[0], right=1.0)
Hof = lambda a: H0*np.sqrt(Om/a**3 + OL)
def lam(a, lam3, lam0):
    """PHYSICAL crossover (kpc). lam3, lam0 are COMOVING crossovers at z = 3 and z = 0. Before z = 3 the comoving crossover
    scales as a^(1/2) (constant adiabatic sound speed, as in our CLASS runs); after z = 3 it grows as a power of a to lam0."""
    if a <= 0.25: lc = lam3*(a/0.25)**0.5
    else: lc = lam3*(lam0/lam3)**(np.log(a/0.25)/np.log(4))
    return lc*a
# ---------------- Part A: linear growth
def growth(kc_Mpc, lam3, lam0):
    kc = kc_Mpc/1e3                                                   # comoving 1/kpc
    def rhs(lna, y):
        a = np.exp(lna); H = Hof(a); dlnH = -1.5*Om/a**3/(H/H0)**2
        mu = 2*np.pi/lam(a, lam3, lam0); kp = kc/a; Sf = mu**2/(kp**2 + mu**2) if np.isfinite(mu) else 1.0
        df, dfp, db, dbp = y; c = 1.5*(H0/H)**2/a**3
        return [dfp, -(2 + dlnH)*dfp + c*Sf*(Oc*df + Ob*db), dbp, -(2 + dlnH)*dbp + c*(Oc*df + Ob*db)]
    a0_ = 1e-3; s = solve_ivp(rhs, [np.log(a0_), 0], [a0_]*4, rtol=1e-7, dense_output=True)
    return s
def ratios(lam3, lam0):
    def at(kc, z, l3, l0):
        s = growth(kc, l3, l0); y = s.sol(np.log(1/(1+z))); return (Oc*y[0] + Ob*y[2])/Om
    ks = np.logspace(-2.5, 1, 40); q = ks/(Om*h**2/h)            # k in 1/Mpc; BBKS with Gamma = Om h (k in h/Mpc -> /h)
    qq = ks/h/(Om*h); T = np.log(1+2.34*qq)/(2.34*qq)*(1+3.89*qq+(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**-0.25
    P = ks**0.965*T**2; W = lambda x: 3*(np.sin(x)-x*np.cos(x))/x**3; R8 = 8/h
    g = np.array([at(kk, 0, lam3, lam0) for kk in ks]); g0 = np.array([at(kk, 0, 1e-9, 1e-9) for kk in ks])
    s8 = np.sqrt(np.trapezoid(P*(g*W(ks*R8))**2*ks**2, ks)/np.trapezoid(P*(g0*W(ks*R8))**2*ks**2, ks))
    lya = (at(5.0, 3, lam3, lam0)/at(5.0, 3, 1e-9, 1e-9))**2
    return s8, lya
# ---------------- Part B: one galaxy
def tz(z): return quad(lambda a: 1/(a*Hof(a)), 0, 1/(1+z))[0]
_ag = np.linspace(2e-3, 1.0, 600); _tg = np.array([tz(1/a - 1) for a in _ag]); aa = lambda t: np.interp(t, _tg, _ag)
hx_ = (1 - fs)*xs; ipk_ = np.argmax(hx_); hmax_ = hx_[ipk_]
def hinv(S):                                                         # solve (1 - f(g/a0)) g = S on the small-g branch (clamped)
    x = np.interp(np.minimum(np.abs(S)/a0, hmax_), hx_[:ipk_+1], xs[:ipk_+1]); return np.sign(S)*x*a0
def solve_Xi(rs_sorted, dm_shell, mu, Xi_guess=None, it=12):
    """Returns Xi' at the shells. Integrated form: r^2 w Xi' = G dM(<r) + mu^2 I(r), I = int r^2 Xi dr, w = 1 - f(|Xi'|/a0).
    Xi' is taken from this integrated form (no numerical differentiation of Xi, which blew up where shells pile up);
    Xi itself (needed only for I) comes from a banded finite-volume solve with w frozen, iterated."""
    r = rs_sorted; N = len(r); re = np.r_[0, 0.5*(r[1:] + r[:-1]), r[-1]*1.5]
    vol_w = (re[1:]**3 - re[:-1]**3)/3; dMenc = np.cumsum(dm_shell)
    Xi = np.zeros(N) if Xi_guess is None else Xi_guess.copy()
    for _ in range(it):
        I = np.cumsum(Xi*vol_w); dXi = hinv((G*dMenc + mu**2*I)/r**2)
        w = np.maximum(1 - f_of(np.abs(dXi)/a0), 1e-4)
        wf = np.r_[0, 0.5*(w[1:] + w[:-1]), w[-1]]
        dr = np.maximum(np.diff(r), 1e-3)
        cpl = np.r_[re[1:-1]**2*wf[1:-1]/dr, 0]
        ab = np.zeros((3, N)); ab[1] = -(cpl + np.r_[0, cpl[:-1]]) - mu**2*vol_w
        ab[0, 1:] = cpl[:-1]; ab[2, :-1] = cpl[:-1]; ab[1, -1] = 1; ab[2, -2] = 0
        rhs = G*dm_shell.copy(); rhs[-1] = 0
        Xi = 0.5*Xi + 0.5*solve_banded((1, 1), ab, rhs)
    I = np.cumsum(Xi*vol_w); return Xi, hinv((G*dMenc + mu**2*I)/r**2)
def galaxy(lam3, lam0, Mb=1e11, N=300, zi=30.0, Mhalo=0.0, kappa=0.0, env=False):
    ai = 1/(1+zi); R = np.geomspace(20, 8000, N); edges = np.r_[0, R]
    m_f = Oc*rhoc*4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
    # optional pre-existing overdensity: a region of total mass Mhalo that would collapse at z ~ 2 in LCDM
    RL = (3*Mhalo/(4*np.pi*Om*rhoc))**(1/3) if Mhalo else 1.0
    Delta = (1.686*(1/31)/(1/3))/(1 + (R/RL)**3) if Mhalo else 0*R
    if env and Mhalo:                                                   # correlated surroundings: mean enclosed overdensity ~ R^-1.2
        Delta = (1.686*(1/31)/(1/3))*np.minimum(1.0, (R/RL)**-1.2)
    r = ai*R*(1 - Delta/3); v = Hof(ai)*r*(1 - Delta/3);
    jang = kappa*np.sqrt(G*Om*rhoc*4/3*np.pi*R**3*R)                   
    return run_shells(R, m_f, r, v, lam3, lam0, Mb, zi, jang)
def run_shells(R, m_f, r, v, lam3, lam0, Mb, zi, jang=0.0):
    N = len(R)
    t = tz(zi); eps = 2.0; Xi = None; out = {}   # jang: specific angular momentum (secondary-infall style)
    marks = [(tz(2.0), "z=2"), (tz(0.9), "z=0.9"), (tz(0.0), "z=0")]
    while t < marks[-1][0]:
        a = aa(t); mu = 2*np.pi/lam(a, lam3, lam0)
        o = np.argsort(r); rs = r[o] + 1e-6*np.arange(N); Mf = np.cumsum(m_f[o])
        rb = np.r_[0, 0.5*(rs[1:] + rs[:-1]), rs[-1]*1.5]
        dm = m_f[o] - Oc*rhoc/a**3*4/3*np.pi*(rb[1:]**3 - rb[:-1]**3)
        if NOXI: dXi = np.zeros(N)
        else: Xi, dXi = solve_Xi(rs, dm, mu, Xi, it=6)
        rr = np.sqrt(rs**2 + eps**2)
        acc_s = -G*(Mb + Mf)/rr**2 + OL*H0**2*rs + dXi
        acc = np.empty(N); acc[o] = acc_s; acc += jang**2/np.maximum(r, 2.0)**3
        dt = max(min(0.01, 0.02*np.min(rr/(np.abs(v[o]) + 1e-3))), DTMIN)
        v += acc*dt; r += v*dt; t += dt
        lo = r < 2.0; r[lo] = 2.0; v[lo] = np.abs(v[lo])
        for tm, lab in marks:
            if t >= tm and lab not in out:
                o = np.argsort(r); out[lab] = (r[o].copy(), np.cumsum(m_f[o]), aa(t))
    return out
if __name__ == "__main__":
    cases = [(100.0, 100.0), (100.0, 1000.0), (100.0, 3000.0), (200.0, 3000.0), (100.0, 10000.0)]
    import os
    if os.environ.get("ONLY_C"):
        cases = [(100.0, 100.0), (100.0, 1000.0), (100.0, 3000.0), (100.0, 10000.0)]
        print("Part C: galaxy (Mb = 1e11) in a pre-existing 1.5e12 Msun overdensity; fluid excess at 30/100/300 kpc (MOND 2.3e11 / 8.8e11 / 2.7e12)")
        for l3, l0 in cases:
            res = galaxy(l3, l0, Mhalo=1.5e12); row = []
            for lab in ("z=2", "z=0.9", "z=0"):
                rs, Mf, av = res[lab]; row.append([np.interp(R_, rs, Mf) - 4/3*np.pi*R_**3*Oc*rhoc/av**3 for R_ in (30, 100, 300)])
            print(f"   lambda {l3:5.0f} -> {l0:5.0f}: " + " | ".join(f"{lab} " + " ".join(f"{x:8.1e}" for x in rr)
                  for lab, rr in zip(("z=2", "z=0.9", "z=0"), row))); sys.stdout.flush()
        raise SystemExit
    print("Part A (linear): comoving crossover (kpc) at z=3 -> today | sigma8/CDM | Lyman-alpha power z=3, k=5/Mpc /CDM")
    for l3, l0 in cases:
        s8, ly = ratios(l3, l0); print(f"   {l3:6.0f} -> {l0:6.0f} | {s8:.3f} | {ly:.3f}"); sys.stdout.flush()
    print("   (allowed: sigma8 within ~5-8% below; Lyman-alpha >= ~0.98)")
    print("\nPart C: same, but the galaxy sits in a pre-existing overdensity of 1.5e12 Msun (an LCDM-like halo region, collapsing ~z 2)")
    for l3, l0 in cases:
        res = galaxy(l3, l0, Mhalo=1.5e12); row = []
        for lab in ("z=2", "z=0.9", "z=0"):
            rs, Mf, av = res[lab]
            row.append([np.interp(R_, rs, Mf) - 4/3*np.pi*R_**3*Oc*rhoc/av**3 for R_ in (30, 100, 300)])
        print(f"   lambda {l3:5.0f} -> {l0:5.0f}: " + " | ".join(f"{lab} " + " ".join(f"{x:8.1e}" for x in rr)
              for lab, rr in zip(("z=2", "z=0.9", "z=0"), row))); sys.stdout.flush()
    print("\nPart B (Mb = 1e11): fluid excess vs MOND halo (MOND: 30 kpc 2.3e11, 100 kpc 8.8e11, 300 kpc 2.7e12)")
    for l3, l0 in cases:
        res = galaxy(l3, l0); row = []
        for lab in ("z=2", "z=0.9", "z=0"):
            rs, Mf, av = res[lab]
            row.append([np.interp(R_, rs, Mf) - 4/3*np.pi*R_**3*Oc*rhoc/av**3 for R_ in (30, 100, 300)])
        print(f"   lambda {l3:5.0f} -> {l0:5.0f}: " + " | ".join(f"{lab} " + " ".join(f"{x:8.1e}" for x in rr)
              for lab, rr in zip(("z=2", "z=0.9", "z=0"), row))); sys.stdout.flush()
