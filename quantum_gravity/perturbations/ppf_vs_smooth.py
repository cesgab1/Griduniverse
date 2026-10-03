"""
Does the treatment of dark-energy FLUCTUATIONS matter for Claim 1?
Our CAMB fits used PPF (an ad hoc smooth recipe that lets w cross -1). In the ghost-free home of the law (VCDM /
extended cuscuton) dark energy has no fluctuations of its own. Proxy: CAMB PPF with its dark-energy perturbations
switched off (internal no_perturbations flag; background only), same background. Compare CMB TT/EE/lensing spectra and sigma8 with PPF.
Metric: cosmic-variance chi2 = sum_l f_sky (2l+1)/2 (dC_l/C_l)^2, f_sky = 0.7 (an upper bound on what any survey can see).
"""
import numpy as np, camb
from scipy.integrate import solve_ivp
H0, ombh2, omch2, mnu = 67.9, 0.02237, 0.1200, 0.06
def law_table(beta=0.5):
    h = H0/100; om = (ombh2 + omch2 + mnu/93.14)/h**2; orad = 4.18e-5/h**2; ol = 1 - om - orad
    g = np.linspace(0.0, np.log(1e-5), 800)
    rhs = lambda lna, y: [beta*(om*np.exp(-3*lna)/2 + orad*np.exp(-4*lna) - np.exp(y[0]))/(om*np.exp(-3*lna) + orad*np.exp(-4*lna) + np.exp(y[0]))]
    s = solve_ivp(rhs, [0, g[-1]], [np.log(ol)], t_eval=g, rtol=1e-9, atol=1e-11)
    a = np.exp(s.t); rde = np.exp(s.y[0]); rm, rr = om*a**-3, orad*a**-4
    return a[::-1], (-1 - beta*(rm/2 + rr - rde)/(rm + rr + rde)/3)[::-1]
a, w = law_table()
def run(model):
    p = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=2600, lens_potential_accuracy=1)
    p.WantTransfer = True; p.set_matter_power(redshifts=[0.0], kmax=2.0)
    if model == "ppf": p.set_dark_energy_w_a(a, w, dark_energy_model="ppf")
    elif model == "smooth":
        p.set_dark_energy_w_a(a, w, dark_energy_model="ppf"); setattr(p.DarkEnergy, "__no_perturbations", True)   # PPF background, DE fluctuations off
        assert getattr(p.DarkEnergy, "__no_perturbations")
    else: pass   # LCDM
    r = camb.get_results(p); cl = r.get_cmb_power_spectra(p, CMB_unit="muK", raw_cl=True)
    return cl["total"], cl["lens_potential"][:, 0], r.get_sigma8_0(), r
out = []
try:
    tP, lP, sP, _ = run("ppf"); tS, lS, sS, _ = run("smooth"); tL, lL, sL, _ = run("lcdm")
except Exception as e:
    out.append(f"CAMB error: {e}"); raise
l = np.arange(len(tP)); fs = 0.7
def cv(x, y, col, lo=2, hi=2500):
    m = (l >= lo) & (l <= hi); return float(np.sum(fs*(2*l[m] + 1)/2*((x[m, col] - y[m, col])/y[m, col])**2))
def cvp(x, y, lo=2, hi=2000):
    m = (l >= lo) & (l <= hi); return float(np.sum(fs*(2*l[m] + 1)/2*((x[m] - y[m])/y[m])**2))
out.append(f"background: Claim 1 law (beta 1/2), H0 {H0}, ombh2 {ombh2}, omch2 {omch2}")
out.append(f"sigma8: PPF {sP:.5f}   no DE fluctuations {sS:.5f}   (diff {100*(sS/sP - 1):+.3f}%)   LCDM same early universe {sL:.5f}")
for lab, col in (("TT", 0), ("EE", 1)):
    out.append(f"{lab}: cosmic-variance chi2 (PPF vs no-fluct) = {cv(tS, tP, col):.3f}   [low l 2-30 only: {cv(tS, tP, col, 2, 30):.3f}]   "
               f"for scale: law vs LCDM = {cv(tP, tL, col):.1f}")
out.append(f"lensing phi-phi: cosmic-variance chi2 (PPF vs no-fluct) = {cvp(lS, lP):.3f}   law vs LCDM = {cvp(lP, lL):.1f}")
m = (l >= 2) & (l <= 30); out.append(f"largest TT difference PPF vs no-fluct (l <= 30): {100*np.max(np.abs(tS[m, 0]/tP[m, 0] - 1)):.3f}%")
mL = (l >= 8) & (l <= 400)
out.append(f"lensing phi-phi difference over the measured range L 8-400 (Planck/ACT lensing): max {100*np.max(np.abs(lS[mL]/lP[mL] - 1)):.2f}%, mean {100*np.mean(lS[mL]/lP[mL] - 1):+.2f}%;  measured lensing amplitude precision ~2-2.5%")
out.append("Caveat: switching DE fluctuations off in CAMB's gauge is a proxy for VCDM (whose momentum constraint differs), not VCDM itself.")
out.append("Reading: chi2 << 1 means no experiment, even a perfect one, could tell the two treatments apart.")
txt = "\n".join(out); print(txt); open("ppf_vs_smooth.txt", "w").write(txt + "\n")
