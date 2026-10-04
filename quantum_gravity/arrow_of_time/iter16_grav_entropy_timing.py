"""
ITERATION 16 (Penrose / Weyl curvature hypothesis): does dark energy's turnover track the production of GRAVITATIONAL
entropy? Several standard measures are computed; NONE is chosen to fit. Background: flat LCDM, Om = 0.31 (the timing of
clumping is set by the background, nearly the same for our law).
 M1  Weyl/Ricci ratio (Wainwright-Goode style), linear theory: tidal (electric Weyl) amplitude E ~ 4 pi G rho_m delta,
     Ricci ~ rho_tot  ->  P(a) ∝ Om(a) D(a). Production = dP/dt.
 M2  Mass fraction in collapsed halos (Press-Schechter), halos of 1e12 Msun and 1e14 Msun: dF/dt.
 M3  Black-hole entropy (Bekenstein-Hawking, the largest gravitational entropy, Penrose): S ∝ sum M_BH^2. BH mass density
     follows the measured cosmic star-formation history (Madau & Dickinson 2014 fit, BH accretion ∝ SFR); if BH masses
     grow in proportion, S ∝ rho_BH^2 -> dS/dt ∝ rho_BH * accretion rate.
 M4  Cosmic-horizon (apparent horizon) entropy S_H ∝ 1/H^2: dS_H/dt.
Compare peak redshifts of the production rates with z = 0.68 (Claim 1 crossing = acceleration onset) and z = 0.46 (toy).
"""
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.special import erfc
Om, OL = 0.31, 0.69
a = np.geomspace(1/30, 1.0, 3000); z = 1/a - 1; H = np.sqrt(Om*a**-3 + OL)
t = cumulative_trapezoid(1/(a*H), a, initial=0)
D = H*cumulative_trapezoid(1/(a*H)**3, a, initial=0) + H*quad(lambda x: 1/(x*np.sqrt(Om*x**-3 + OL))**3, 1e-8, a[0])[0]
D /= D[-1]; Oma = Om*a**-3/H**2
def peak(rate): i = np.argmax(rate); return z[i]
out = ["ITERATION 16: when is gravitational entropy produced fastest?  (peak redshift of each production rate)", ""]
P = Oma*D; r1 = np.gradient(P, t); out.append(f"M1 Weyl/Ricci (linear): P ∝ Om(a) D(a): P peaks at z = {peak(P):.2f}; production dP/dt peaks at z = {peak(r1):.2f}, turns negative at z = {z[np.where(np.diff(np.sign(r1)))[0][0]] if np.any(np.diff(np.sign(r1))) else float('nan'):.2f}")
for M, sig0 in ((1e12, 2.0), (1e14, 0.8)):              # sigma(M) today for these masses (approx., sigma8 = 0.81)
    nu = 1.686/(sig0*D); F = erfc(nu/np.sqrt(2)); r2 = np.gradient(F, t)
    out.append(f"M2 halo collapse, M = {M:.0e} Msun: dF/dt peaks at z = {peak(r2):.2f}")
sfr = 0.015*(1 + z)**2.7/(1 + ((1 + z)/2.9)**5.6)       # Msun/yr/Mpc^3, Madau & Dickinson 2014
dt_yr = np.gradient(t)*14.5e9/1.0                        # t in Hubble times -> yr (h = 0.674)
rho_bh = cumulative_trapezoid(sfr*np.gradient(t)*14.5e9, initial=0)  # ∝ BH mass density (BH accretion ∝ SFR)
r3 = rho_bh*sfr; out.append(f"M3 black-hole entropy: accretion (SFR) peaks at z = {peak(sfr):.2f}; entropy production rho_BH x rate peaks at z = {peak(r3):.2f}")
SH = 1/H**2; r4 = np.gradient(SH, t); out.append(f"M4 horizon entropy 1/H^2: production peaks at z = {peak(r4):.2f} (still rising today: {r4[-1] > r4[-2]})")
out += ["", "Reference: acceleration onset (q = 0, Claim 1 crossing) z = 0.68; toy (memory kappa = 3) crossing z = 0.46.",
        "VERDICT: the measures disagree (peaks from z ~ 0 to 1.4). The only ones near q = 0 (M1 turnover 0.79, M4 0.65) are",
        "driven by the expansion history itself (M4: d(1/H^2)/dt = 2(1+q)/H is a function of q), so the match is circular,",
        "not independent evidence. NO physical link between dark energy's turnover and gravitational-entropy production found."]
txt = "\n".join(out); print(txt); open("iter16_grav_entropy_timing.txt", "w").write(txt + "\n")
