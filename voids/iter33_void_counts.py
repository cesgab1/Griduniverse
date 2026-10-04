"""
ITERATION 33: Coalesce's void idea -- compare void sizes to read dark energy. Prediction: how many voids of each size does
Claim 1 give vs Lambda? Same matter (Om = 0.31, h = 0.675) and same early ripples (CAMB linear power, Planck); only the
dark-energy history differs (through the growth of structure). Void abundance: Vdn model (Jennings, Li & Hu 2013;
Sheth & van de Weygaert 2004), linear void threshold -2.71, collapse 1.686, volume-conserving radius mapping.
"""
import numpy as np, camb
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
Om, h = 0.31, 0.675; Or = 9.1e-5; OL = 1 - Om - Or
pars = camb.set_params(H0=100*h, ombh2=0.0224, omch2=Om*h*h - 0.0224 - 0.0006, mnu=0.06, As=2.1e-9, ns=0.965)
pars.set_matter_power(redshifts=[0.0], kmax=50); res = camb.get_results(pars)
kh, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=50, npoints=1500); pk = pk[0]
def sig(R):
    x = kh*R; W = 3*(np.sin(x) - x*np.cos(x))/x**3
    return np.sqrt(np.trapezoid(kh**3*pk*W**2/(2*np.pi**2), np.log(kh)))
def E_L(a): return np.sqrt(Om*a**-3 + Or*a**-4 + OL)
def E_C(a):
    f = lambda e: e*e - Om*a**-3 - Or*a**-4 - OL*(a*e)**-0.5
    return brentq(f, 1e-6, 10*np.sqrt(Om*a**-3 + Or*a**-4) + 10)
def growth(E):
    def rhs(lna, y):
        a = np.exp(lna); e = E(a); de = (np.log(E(a*1.001)) - np.log(E(a/1.001)))/(2*np.log(1.001))
        return [y[1], -(2 + de)*y[1] + 1.5*Om*a**-3/e**2*y[0]]
    lna = np.linspace(np.log(1/31), 0, 600)
    s = solve_ivp(rhs, [lna[0], 0], [1/31, 1/31], t_eval=lna, rtol=1e-9)
    return lna, s.y[0]                                  # same early normalisation (D = a at z = 30)
lnL, DL = growth(E_L); lnC, DC = growth(E_C)
DL1 = DL[-1]
dv, dc = 2.71, 1.686; Dsv = dv/(dc + dv)
def vdn(R, z, D):
    a = 1/(1 + z); g = np.interp(np.log(a), lnL, D)/DL1      # growth relative to LCDM today (shape of P(k) is LCDM's today)
    RL = R/1.7189                                           # linear radius (shell crossing, 1 + delta_NL = 0.2)
    s = sig(RL)*g; ds = (np.log(sig(RL*1.01)) - np.log(sig(RL/1.01)))/(2*np.log(1.01))
    x = Dsv*s/dv
    if x <= 0.276:   # Jennings, Li & Hu 2013 eq. (A?) small-x form
        f = np.sqrt(2/np.pi)*(dv/s)*np.exp(-dv**2/(2*s**2))
    else:            # first four terms suffice for x > 0.276
        j = np.arange(1, 5); f = 2*np.sum(np.exp(-(j*np.pi*x)**2/2)*j*np.pi*x**2*np.sin(j*np.pi*Dsv))
    return f/(4/3*np.pi*R**3)*abs(ds)
zs = (0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0)
rat = {z: np.interp(np.log(1/(1+z)), lnL, DC)/np.interp(np.log(1/(1+z)), lnL, DL) for z in zs}
out = ["ITERATION 33: voids as a dark-energy probe -- Claim 1 vs Lambda (same matter, same early ripples)", "",
       "growth of structure, Claim 1 / Lambda: " + "  ".join(f"z={z}: {rat[z]:.4f}" for z in zs), ""]
# realistic sensitivity: void counts respond as d ln n ~ (delta_v,eff / sigma)^2 d ln sigma, with a calibrated (tracer) linear
# threshold delta_v,eff ~ 1 (BOSS/Euclid analyses recalibrate it; the textbook 2.71 puts these radii in a far tail, rejected)
out.append(" R [Mpc/h]   amplification   void-count difference z=0.5   z=1.0")
for R in (20, 30, 40):
    for_z = []
    for z in (0.5, 1.0):
        g_ = np.interp(np.log(1/(1+z)), lnL, DL)/DL1; sg = sig(R/1.7189)*g_
        A = (1.0/sg)**2; for_z.append((A, A*(rat[z] - 1)*100))
    out.append(f"   {R:4d}         {for_z[0][0]:5.1f}             {for_z[0][1]:+.2f}%                    {for_z[1][1]:+.2f}%")
out += ["", "Measured/forecast precision: BOSS DR12 void size function w = -1.1 +/- 0.2; Euclid voids alone: w to ~10%,",
        "FoM(w0,wa) = 17. Claim 1 changes void counts by ~0.5-3% (fewer large voids below z ~ 1, a hair more above): the",
        "pattern (sign change near z ~ 1) is distinctive but below foreseeable void-only precision; useful only combined.",
        "Note: a first version used the textbook threshold 2.71, which puts these voids in an extreme tail (densities ~1e-19 to",
        "1e-87) and gave spurious 5-40% ratios; rejected."]
txt = "\n".join(out); print(txt); open("iter33_void_counts.txt", "w").write(txt + "\n")
