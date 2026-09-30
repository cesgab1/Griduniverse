"""
One mix for clusters: after the grid switch (z_s, derived: gas decouples from the CMB), a fraction f of the dark ingredient keeps
gathering into structures like cold dark matter; the rest (1-f) stays smooth (part of the grid, no clumping). Slack acts on bound matter.
Inputs derived or taken from our own fits (no hand-set numbers):
  cosmology  : forward CMB fit (fwd_gridE043.json: H0, ombh2, omch2, As, ns)
  a0         : c H0 / 6 (elastic grid)
  z_s        : 138 (gas-CMB halfway decoupling, gas_decoupling.py)
  cluster dark/baryon ratio : cosmic ratio f*Ω_x/Ω_b (clusters as fair samples)
Data (external, measured): cluster abundance consistent with Planck-ΛCDM (eROSITA eRASS1 counts, S8 = 0.86 ± 0.01);
  per-cluster pull: CLASH/BCG radial acceleration relation log g_obs = 0.52 log g_bar - 4.19 (Tian+2024).
"""
import json, numpy as np, camb
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import erfc
fit = json.load(open("/home/claude/fullfit/pact/fwd_gridE043.json"))["params"]
G=6.674e-11; Msun=1.989e30; Mpc=3.0857e22; c=2.998e8
h=fit["H0"]/100; H0=fit["H0"]*1e3/Mpc; wb, wc = fit["ombh2"], fit["omch2"]
Ob, Ox = wb/h**2, wc/h**2; Om = Ob+Ox; OL = 1-Om
a0 = c*H0/6; zs = 138.0; as_ = 1/(1+zs)
rhoc = 3*H0**2/(8*np.pi*G); rho_m0 = Om*rhoc
nu = lambda y: 1/(-np.expm1(-np.sqrt(y)))
print(f"derived: H0 {fit['H0']:.2f}, Ω_b {Ob:.4f}, Ω_x {Ox:.4f}, a0 = cH0/6 = {a0:.3e}, switch z = {zs:.0f}")
def t_of_a(a): return 2/(3*H0*np.sqrt(OL))*np.arcsinh(np.sqrt(OL/Om)*a**1.5)
def a_of_t(t): return (np.sinh(1.5*H0*np.sqrt(OL)*t)*np.sqrt(Om/OL))**(2/3)
LAM = 3*H0**2*OL
def collapse_a(M, di, f, slack, a_i=1e-3):
    """shell enclosing total matter M (Lagrangian). Before a_s all matter clusters; after, the unclustered dark part (1-f)Ω_x/Ω_m of the
    enclosed dark mass is replaced by a smooth background share."""
    r_i = (3*M/(4*np.pi*rho_m0*a_i**-3*(1+di)))**(1/3); t_i = t_of_a(a_i); Hi = H0*np.sqrt(Om*a_i**-3+OL)
    v_i = Hi*r_i*(1-di/3); t_s = t_of_a(as_); smooth = (1-f)*Ox/Om
    def Menc(t, r):
        if t < t_s: return M
        rho_bg = rho_m0*a_of_t(t)**-3
        return (1-smooth)*M + smooth*rho_bg*4*np.pi/3*r**3
    def rhs(t, y):
        r, v = y; g = G*Menc(t, r)/r**2
        if slack and v < 0: g *= nu(g/a0)
        return [v, -g + LAM/3*r]
    tend = t_of_a(2.0)
    e1 = lambda t, y: y[1]; e1.direction = -1
    s = solve_ivp(rhs, [t_i, tend], [r_i, v_i], events=e1, rtol=1e-8, atol=1e-6*r_i, dense_output=True, max_step=tend/3000)
    if not len(s.t_events[0]): return 10
    tta = s.t_events[0][0]; rta = s.sol(tta)[0]
    e2 = lambda t, y: y[0]-rta/2; e2.terminal = True; e2.direction = -1
    s2 = solve_ivp(rhs, [tta, tend], s.sol(tta), events=e2, rtol=1e-8, atol=1e-6*r_i, max_step=(tend-tta)/4000)
    return a_of_t(s2.t_events[0][0]) if len(s2.t_events[0]) else 10
def thr(M, f, slack): return brentq(lambda d: collapse_a(M, d, f, slack) - 1.0, 1e-5, 0.2, xtol=1e-8)
p = camb.set_params(H0=fit["H0"], ombh2=wb, omch2=wc, As=np.exp(fit["logA"])*1e-10, ns=fit["ns"], mnu=0.06)
p.set_matter_power(redshifts=[0], kmax=50); r = camb.get_results(p)
Ms = np.array([1e14, 3e14, 1e15])
Rs = (3*Ms*Msun/(4*np.pi*rho_m0))**(1/3)/Mpc*h
sig = r.get_sigmaR(Rs, hubble_units=True)[0] if np.ndim(r.get_sigmaR(Rs, hubble_units=True)) > 1 else r.get_sigmaR(Rs, hubble_units=True)
print(f"σ(M) today: " + ", ".join(f"{M:.0e}: {s:.3f}" for M, s in zip(Ms, sig)), f"  σ8 = {r.get_sigma8_0():.3f}")
base = [thr(M*Msun, 1.0, False) for M in Ms]          # ΛCDM (f=1, no slack)
# per-cluster pull: CLASH relation vs model g = gN ν(gN/a0), gN = g_bar (1 + f Ω_x/Ω_b)
gbar = np.logspace(-12, -9, 13); gmeas = 10**(0.52*np.log10(gbar) - 4.19)
def pull_chi(f):
    gN = gbar*(1 + f*Ox/Ob); gm = gN*nu(gN/a0)
    return np.sum(((np.log10(gm) - np.log10(gmeas))/0.05)**2), np.median(gm/gmeas)
out = []
print(f"\n{'f':>5} | cluster counts vs observed (≈ΛCDM), M = 1e14, 3e14, 1e15 Msun | per-cluster pull model/measured (median)  χ²(13 pts, 0.05 dex)")
for f in [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0]:
    ratios = []
    for M, b, s in zip(Ms, base, sig):
        d = thr(M*Msun, f, True); dc = 1.686*d/b
        ratios.append(erfc(dc/np.sqrt(2)/s)/erfc(1.686/np.sqrt(2)/s))
    chi, med = pull_chi(f)
    out.append({"f": f, "counts": ratios, "pull_median": med, "pull_chi2": chi})
    print(f"{f:5.2f} | " + "  ".join(f"×{x:7.3f}" for x in ratios) + f"   | ×{med:5.2f}   χ² {chi:8.1f}")
json.dump(out, open("cluster_mix.json", "w"), indent=1)
