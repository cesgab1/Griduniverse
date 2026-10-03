"""
ITERATION 8: how faint must the leftover dark-energy lumps be? Compare their large-angle CMB imprint with the real sky.

Lumps (iterations 5-7): delta = relative dark-energy fluctuation, Gaussian for many layers (central limit), with
  spatial correlation   <delta delta'> = (2/N_eff) f(d/L),  f(x) = ((1 - e^-x)/x)^2   (quadratic energy: square of the
                        tension correlation (1 - e^-x)/x), comoving L = c/(kappa a H), kappa = 3
  time correlation      exp(-lambda |t - t'|), lambda = 2 kappa H  (energy decorrelates twice as fast as the tension)
Potential (comoving Poisson, units c = H0 = 1):  Phi_k = (3/2) Omega_DE a^2 delta_k / k^2   (sub-horizon form; matter's
response to the lumps neglected: both are order-1 caveats, fine for a bound on N within a factor of a few)
CMB imprint: late ISW, Theta = 2 ∫ dt dPhi/dt along the line of sight = -2 ∫ dt [d j_l(k chi)/dt] Phi  (boundary terms vanish)
  C_l = (2/pi) ∫ k^2 dk  4 ∫∫ dt dt' J(t) J(t') A(t) A(t') P(k; t, t') e^{-lambda |t - t'|}
Real sky: Planck-like LCDM spectrum from CAMB; cosmic variance sigma_l = C_l sqrt(2/(2l+1)) for l = 2..30.
Bound: added power must not exceed chi^2 = 4 (about 2 sigma) summed over l = 2..30.
"""
import numpy as np, camb
from scipy.integrate import quad, cumulative_trapezoid
from scipy.special import spherical_jn
Om, OD, kappa = 0.31, 0.69, 3.0
f = lambda x: ((1 - np.exp(-x))/x)**2 if x > 1e-8 else 1.0
# F(q) = ∫ f(x) x^2 sinc(q x) dx = (1/q) ∫ f(x) x sin(q x) dx   (oscillatory, Fourier-weight quadrature)
qs = np.geomspace(1e-3, 1e3, 241)
Fq = np.array([quad(lambda x: f(x)*x, 0, np.inf, weight="sin", wvar=q, limlst=200)[0]/q for q in qs])
F = lambda q: np.exp(np.interp(np.log(q), np.log(qs), np.log(np.abs(Fq))))
# background (flat LCDM approximation for the lump transport)
a = np.geomspace(0.08, 1.0, 260); H = np.sqrt(Om*a**-3 + OD)
t = cumulative_trapezoid(1/(a*H), a, initial=0)                     # time since a = 0.08 (units 1/H0)
chi = np.array([quad(lambda x: 1/(x*x*np.sqrt(Om*x**-3 + OD)), ai, 1)[0] for ai in a])
Lc = 1/(kappa*a*H); lam = 2*kappa*H; A_a = 1.5*OD*a**2           # A(t)/k^2 has the 1/k^2 put in below
dt = np.gradient(t)
def Cl_extra(ell, Neff=1.0):
    ks = np.geomspace(0.05, 60, 220); tot = np.zeros_like(ks)
    for i, k in enumerate(ks):
        jl = spherical_jn(ell, k*chi); J = np.gradient(jl, t)          # d j_l(k chi(t))/dt
        g = J*A_a/k**2*np.sqrt(dt)
        Lm = np.sqrt(np.outer(Lc, Lc)); lm = 0.5*(lam[:, None] + lam[None, :])
        P = 4*np.pi*Lm**3*(2/Neff)*F(np.clip(k*Lm, 1e-3, 1e3))
        Kt = np.exp(-lm*np.abs(t[:, None] - t[None, :]))
        tot[i] = 4*np.einsum("i,ij,j->", g, P*Kt, g)
    return 2/np.pi*np.trapezoid(ks**2*tot, ks)
pars = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.120, mnu=0.06, tau=0.054, As=2.1e-9, ns=0.965, lmax=60)
res = camb.get_results(pars); cl = res.get_cmb_power_spectra(pars, CMB_unit=None, raw_cl=True)["total"][:, 0]   # dimensionless (dT/T)^2
ells = np.arange(2, 31); ce = np.array([Cl_extra(l) for l in ells]); base = cl[ells]
sig = base*np.sqrt(2/(2*ells + 1)); chi2_1 = np.sum((ce/sig)**2)
Nmin = np.sqrt(chi2_1/4.0)
T0 = 2.7255e6
out = ["ITERATION 8: large-angle CMB imprint of the leftover dark-energy lumps (late ISW)", "",
       " l    D_l real sky [muK^2]   D_l lumps for ONE layer [muK^2]"]
for l, b, c in zip(ells[:8], base[:8], ce[:8]):
    out.append(f" {l:2d}    {l*(l+1)*b/2/np.pi*T0**2:10.1f}            {l*(l+1)*c/2/np.pi*T0**2:10.3e}")
out += ["", f"one layer: added power exceeds cosmic variance by chi^2 = {chi2_1:.2e} (l = 2-30)",
        f"-> independent layers needed: N_eff >= {Nmin:.1e}   (2 sigma; added power ∝ 1/N_eff, chi^2 ∝ 1/N_eff^2)",
        f"   compare the rough estimate of iteration 6: 2.6e10.  With neighbour leakage r per gap, N = N_eff (1 + 2 r^2/(1 - r^2)).",
        f"   window from the species bound + gamma-ray-burst timing: N <~ 1.2e15  -> {'OPEN' if Nmin < 1.2e15 else 'CLOSED'}",
        "",
        "Signature if N_eff sits near the bound: extra large-angle power with a characteristic shape (shown above), on top of",
        "LCDM; note the real sky has LOW power at l = 2 (the quadrupole anomaly), so added power there is disfavoured, not favoured."]
txt = "\n".join(out); print(txt); open("iter8_cmb_lumps.txt", "w").write(txt + "\n")
np.save("iter8_cl_one_layer.npy", np.c_[ells, ce, base])
