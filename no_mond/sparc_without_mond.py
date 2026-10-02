"""
Without MOND, galaxies need dark halos (Newtonian gravity: baryons + fluid halo). Required 'absolute truth' #1: do SPARC
rotation curves (and the tight radial acceleration relation they define) come out right?
Models (Newtonian gravity, stellar M/L 0.5 disk / 0.7 bulge as in our MOND fits):
  A  MOND (baryons only, a0 = 1.2e-10, nu = McGaugh)                         0 free numbers per galaxy
  B  fluid = CDM-like NFW halo, mass from abundance matching (Moster+2013), concentration from Dutton & Maccio 2014   0 free
  C  NFW, halo mass free (concentration from the c-M relation)              1 free per galaxy
  D  cored halo (Burkert), mass and core radius free                        2 free per galaxy
  E  MOND with the disk M/L free (bulge 1.4x disk)                         1 free per galaxy
Scores: total chi2 (5% velocity error floor), BIC (chi2 + k ln N per galaxy), and the scatter of the radial acceleration relation
each model produces (observed scatter is the benchmark).
"""
import numpy as np, glob
from scipy.optimize import minimize, minimize_scalar, brentq
G = 4.30e-6; a0 = 1.2e-10*3.086e19/1e6; rhoc = 136.0
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def shmr(M):  # Moster+2013 z=0
    M1, N, b, g = 10**11.59, 0.0351, 1.376, 0.608; return 2*N*M/((M/M1)**-b + (M/M1)**g)
def Mhalo_from_Mstar(Ms): return np.exp(brentq(lambda lm: np.log(shmr(np.exp(lm))) - np.log(Ms), np.log(1e8), np.log(1e16)))
def conc(M): return 10**(0.905 - 0.101*np.log10(M*0.7/1e12))
def nfw_v2(r, M):
    R200 = (3*M/(4*np.pi*200*rhoc))**(1/3); c = conc(M); x = r*c/R200; m = lambda x: np.log(1 + x) - x/(1 + x)
    return G*M*m(x)/m(c)/r
def bur_v2(r, M, rc):
    R200 = (3*M/(4*np.pi*200*rhoc))**(1/3); f = lambda x: np.log(1 + x**2) + 2*np.log(1 + x) - 2*np.arctan(x)
    return G*M*f(r/rc)/f(R200/rc)/r
gals = []; gal_parts = {}
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb, _, _ = d.T
    ok = (r > 0) & (V > 0); r, V, eV, Vg, Vd, Vb = r[ok], V[ok], eV[ok], Vg[ok], Vd[ok], Vb[ok]
    if len(r) < 5: continue
    vb2 = Vg*abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2
    Ms = (0.5*Vd[-1]**2 + 0.7*Vb[-1]**2)*r[-1]/G                     # stellar mass (outermost enclosed, ~0.9 of total)
    e = np.sqrt(eV**2 + (0.05*V)**2)
    gals.append((fn.split("/")[-1][:-11], r, V, e, vb2, max(Ms, 1e6))); gal_parts[fn.split("/")[-1][:-11]] = (Vg, Vd, Vb)
def chi(v2, V, e): return np.sum(((np.sqrt(np.maximum(v2, 0)) - V)/e)**2)
tot = {k: 0.0 for k in "ABCDE"}; bic = {k: 0.0 for k in "ABCDE"}; npts = 0; resid = {k: [] for k in "ABCDE"}; kpar = {"A": 0, "B": 0, "C": 1, "D": 2, "E": 1}
for name, r, V, e, vb2, Ms in gals:
    n = len(r); npts += n
    gN = np.maximum(vb2, 1e-9)/r
    mods = {}
    mods["A"] = nu(gN/a0)*gN*r
    Mh = Mhalo_from_Mstar(Ms); mods["B"] = vb2 + nfw_v2(r, Mh)
    c1 = minimize_scalar(lambda lm: chi(vb2 + nfw_v2(r, 10**lm), V, e), bounds=(8, 14.5), method="bounded")
    mods["C"] = vb2 + nfw_v2(r, 10**c1.x)
    best = None
    for rc0 in (0.5, 2, 8):
        s = minimize(lambda p: chi(vb2 + bur_v2(r, 10**p[0], 10**p[1]), V, e), [11, np.log10(rc0)], method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-3, maxiter=400))
        if best is None or s.fun < best.fun: best = s
    mods["D"] = vb2 + bur_v2(r, 10**best.x[0], 10**best.x[1])
    # E: MOND with the disk mass-to-light ratio free (1 per galaxy), the standard way MOND fits are done
    Vg_, Vd_, Vb_ = gal_parts[name]
    def mond_v2(ml):
        vb2_ = Vg_*abs(Vg_) + ml*Vd_**2 + 1.4*ml*Vb_**2; gN_ = np.maximum(vb2_, 1e-9)/r; return nu(gN_/a0)*gN_*r
    e1 = minimize_scalar(lambda ml: chi(mond_v2(ml), V, e), bounds=(0.1, 1.5), method="bounded")
    mods["E"] = mond_v2(e1.x)
    for k, v2 in mods.items():
        c = chi(v2, V, e); tot[k] += c; bic[k] += c + kpar[k]*np.log(n)
        go_mod = np.maximum(v2, 1e-9)/r; resid[k].extend(np.log10(go_mod/gN))
obs = []
for name, r, V, e, vb2, Ms in gals:
    gN = np.maximum(vb2, 1e-9)/r; obs.extend(np.log10(V**2/r/gN))
obs = np.array(obs)
print(f"{len(gals)} SPARC galaxies, {npts} points (5% error floor)")
lab = {"E": "MOND, disk M/L free (1/galaxy)", "A": "MOND (0 free)", "B": "CDM-like halo, abundance matched (0 free)", "C": "NFW, mass free (1/galaxy)", "D": "cored halo (2/galaxy)"}
for k in "AEBCD":
    print(f"  {lab[k]:42s} chi2 {tot[k]:9.0f}  chi2/pt {tot[k]/npts:6.2f}  BIC {bic[k]:9.0f}")
# RAR scatter: residual of log g_obs about the mean relation, binned in g_bar
gb_all = []
for name, r, V, e, vb2, Ms in gals: gb_all.extend(np.log10(np.maximum(vb2, 1e-9)/r))
gb_all = np.array(gb_all)
def rar_scatter(y):
    y = np.array(y) + gb_all; bins = np.linspace(gb_all.min(), gb_all.max(), 15); s = []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (gb_all >= a) & (gb_all < b)
        if m.sum() > 20: s.append(np.std(y[m]))
    return np.mean(s)
pass  # print(f"  RAR scatter (dex, binned): observed {rar_scatter(obs - 0):.3f} | " + " | ".join(f"{k} {rar_scatter(np.array(resid[k])):.3f}" for k in "ABCD"))

# F: 'salt' rule without MOND (fluid = 5.4 x retained baryons, spread like hot gas: rho ~ r^-2 out to 143 kpc), 0 free
totF = 0.0
for name, r, V, e, vb2, Ms in gals:
    Vg_, Vd_, Vb_ = gal_parts[name]; Mgas = Vg_[-1]**2*r[-1]/G*1.33
    Mret = Ms + Mgas + Ms                                                # stars + cold gas + hot gas (= M*, our KiDS model)
    vf2 = G*5.4*Mret*np.minimum(r, 143.0)/143.0/r
    totF += chi(vb2 + vf2, V, e)
print(f"  {'salt rule without MOND (0 free)':42s} chi2 {totF:9.0f}  chi2/pt {totF/npts:6.2f}")
