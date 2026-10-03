"""
ITERATION 9: let ordinary + dark matter fall into the dark-energy lumps.
Growth with an extra source (physical time, units H0 = c = 1):
   delta_m'' + 2 H delta_m' - (3/2) Om a^-3 delta_m = (3/2) OD delta_DE        (4 pi G rho_DE delta_DE on the right)
Solved with the exact LCDM Green's function (growing D+, decaying D- = H). Linear algebra on a time grid:
   delta_m = G delta_DE ;   Phi = (3/2)[OD a^2 delta_DE + Om a^-1 delta_m] / k^2   (comoving Poisson)
Lumps as in iteration 8 (spatial ((1-e^-x)/x)^2 on c/(3aH), time decorrelation 6H, variance 2/N_eff).
Outputs: (1) revised CMB bound including the matter response; (2) induced matter power today vs LCDM P(k) on the
largest scales surveys measure.
"""
import numpy as np, camb
from scipy.integrate import quad, cumulative_trapezoid
from scipy.special import spherical_jn
Om, OD, kappa = 0.31, 0.69, 3.0
f = lambda x: ((1 - np.exp(-x))/x)**2 if x > 1e-8 else 1.0
qs = np.geomspace(1e-3, 1e3, 241)
Fq = np.array([quad(lambda x: f(x)*x, 0, np.inf, weight="sin", wvar=q, limlst=200)[0]/q for q in qs])
F = lambda q: np.exp(np.interp(np.log(q), np.log(qs), np.log(np.abs(Fq))))
a = np.geomspace(0.08, 1.0, 240); H = np.sqrt(Om*a**-3 + OD)
t = cumulative_trapezoid(1/(a*H), a, initial=0); dt = np.gradient(t)
chi = np.array([quad(lambda x: 1/(x*x*np.sqrt(Om*x**-3 + OD)), ai, 1)[0] for ai in a])
Lc = 1/(kappa*a*H); lam = 2*kappa*H
# Green's function: D- = H, D+ = H ∫ da/(aH)^3 ; G(t,s) = [D+(t)D-(s) - D-(t)D+(s)] / W(s), W = D+' D- - D-' D+
Ip = cumulative_trapezoid(1/(a*H)**3, a, initial=0) + quad(lambda x: 1/(x*np.sqrt(Om*x**-3 + OD))**3, 1e-6, a[0])[0]
Dp, Dm = H*Ip, H
dDp, dDm = np.gradient(Dp, t), np.gradient(Dm, t)
W = dDp*Dm - dDm*Dp
G = (np.outer(Dp, Dm) - np.outer(Dm, Dp))/W[None, :]
G = np.tril(G, -1)*dt[None, :]*1.5*OD                                   # include source factor and quadrature weight
def Mcov(k, Neff):
    Lm = np.sqrt(np.outer(Lc, Lc)); lm = 0.5*(lam[:, None] + lam[None, :])
    return 4*np.pi*Lm**3*(2/Neff)*F(np.clip(k*Lm, 1e-3, 1e3))*np.exp(-lm*np.abs(t[:, None] - t[None, :]))
def Cl(ell, with_matter, Neff=1.0):
    ks = np.geomspace(0.05, 60, 200); tot = np.zeros_like(ks)
    for i, k in enumerate(ks):
        J = np.gradient(spherical_jn(ell, k*chi), t); g = J*np.sqrt(dt)*np.sqrt(dt)    # ∫ dt weights (both sides)
        T = np.diag(1.5*OD*a**2) + (np.diag(1.5*Om/a) @ G if with_matter else 0)
        T = T/k**2; C = T @ Mcov(k, Neff) @ T.T
        gg = J*dt; gg[0] += spherical_jn(ell, k*chi[0]); tot[i] = 4*gg @ C @ gg   # + early boundary term
    return 2/np.pi*np.trapezoid(ks**2*tot, ks)
pars = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.120, mnu=0.06, tau=0.054, As=2.1e-9, ns=0.965, lmax=60)
pars.WantTransfer = True; pars.set_matter_power(redshifts=[0.0], kmax=1.0)
res = camb.get_results(pars); base = res.get_cmb_power_spectra(pars, CMB_unit=None, raw_cl=True)["total"][:, 0]
ells = np.arange(2, 31); sig = base[ells]*np.sqrt(2/(2*ells + 1))
out = ["ITERATION 9: matter falls into the dark-energy lumps", ""]
for wm in (False, True):
    ce = np.array([Cl(l, wm) for l in ells]); chi2 = np.sum((ce/sig)**2)
    out.append(f"CMB (late ISW) one-layer chi2 = {chi2:.2e}   {'WITH' if wm else 'without'} matter response  ->  N_eff >= {np.sqrt(chi2/4):.1e}")
# induced matter power today
h = 0.674; kh = np.geomspace(1e-3, 0.1, 9)
kh_c, z_c, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=1.0, npoints=400)
Nb = None
out.append("")
out.append("Induced matter power today at the CMB bound, vs LCDM linear P(k):")
for N in (3e7, 1e8):
    out.append(f"  N_eff = {N:.0e}:")
    for khv in kh:
        k = khv*h*2997.92458                                              # h/Mpc -> units of H0/c
        C = G @ Mcov(k, N) @ G.T; Pm_extra = C[-1, -1]*(2997.92458/h)**3    # (Mpc/h)^3
        Plcdm = np.interp(khv, kh_c, pk[0]); out.append(f"     k = {khv:6.3f} h/Mpc:  extra/LCDM = {Pm_extra/Plcdm:.1e}")
out.append("Reading: if extra/LCDM << 1 on the largest surveyed scales (k ~ 0.002-0.02 h/Mpc, where survey errors are ~10-50%),")
out.append("galaxy clustering does not tighten the bound; the CMB large-angle sky stays the strongest test.")
txt = "\n".join(out); print(txt); open("iter9_matter_response.txt", "w").write(txt + "\n")
