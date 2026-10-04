"""ITERATION 63 (pre-registered in PREREG_63.md): uncertainty audit of the 24% excess."""
import numpy as np, camb, ctypes
from camb.baseconfig import camblib
from scipy.optimize import brentq
alpha = 1/137.035999084
camblib.set_vconst.argtypes = [ctypes.c_double]*2; camblib.set_vswitch.argtypes = [ctypes.c_double]*2
def history(me):
    camblib.set_vconst(1.0, me); camblib.set_vswitch(-1.0, 30.0)
    p = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.1200, YHe=0.245, tau=0.054); p.Recomb.use_rosenbrock = False
    r = camb.get_results(p); z = np.logspace(1.0, 3.3, 5000)
    ev = r.get_background_redshift_evolution(z, ["x_e"], format="array"); xe = ev[:, 0]
    H = r.hubble_parameter(z)*1e3/3.0857e22; Tc = 2.7255*(1 + z); fHe = 0.245/(4*(1 - 0.245))
    me_kg = 9.109e-31*me
    CH = 8*6.6524e-29/me**2*7.5657e-16*Tc**4/(3*me_kg*2.998e8)*xe/(1 + fHe + xe)/H     # Thomson cross-section ~ 1/m_e^2
    return z, CH
def zs_of(z, CH, k):
    m = z > 30; zz, yy = z[m], CH[m]; i = np.where(np.diff(np.sign(yy - k)))[0]; return float(np.interp(k, [yy[i[-1]], yy[i[-1]+1]], [zz[i[-1]], zz[i[-1]+1]]))
Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or
def E1(zz, beta=0.5):
    a = 1/(1 + zz); g = lambda le: np.exp(2*le) - Om*a**-3 - Or*a**-4 - OL*(a*np.exp(le))**-beta
    return np.exp(brentq(g, -10, 40))
rho_meas = OL*1.0537e4*0.674**2
def ratio(zs, beta=0.5, ombh2=0.0224, Yp=0.245):
    n_e0 = ombh2*1.878e-29/1.6726e-24*(1 - Yp/2)
    return alpha*511e3*n_e0*(1 + zs)**3*(E1(zs, beta)/(1 + zs))**beta/rho_meas
out = ["ITERATION 63: is the 24% excess solid? (decision rule pre-registered)", ""]
z1, C1 = history(1.0); za, Ca = history(1 + alpha)
base = zs_of(z1, C1, 1.0); r0 = ratio(base)
out.append(f"baseline (standard electron, k = 1): payday z = {base:.1f}, prediction/measured = {r0:.3f}")
out.append(f"U1 sensitivity: d ln(ratio) / d ln(1+z_s) = {np.log(ratio(base*1.01+0.01)/r0)/np.log(1.01):.2f}")
for k in (0.5, 1.0, 2.0):
    zk = zs_of(z1, C1, k); out.append(f"U2 criterion k = {k}: payday z = {zk:.1f}, ratio = {ratio(zk):.2f}")
zA = zs_of(za, Ca, 1.0); out.append(f"U3 heavier electron (1 + alpha) in recombination: payday z = {zA:.1f}, ratio = {ratio(zA):.2f}")
for b in (0.41, 0.5, 0.63, 0.85):
    out.append(f"U4 Claim 1 exponent beta = {b}: ratio = {ratio(base, beta=b):.2f}")
out.append(f"U5 omega_b h^2 +/- 0.0001: {ratio(base, ombh2=0.0225)/r0 - 1:+.3f};  Y_p +/- 0.003: {ratio(base, Yp=0.248)/r0 - 1:+.3f};  "
           f"measured omega_DE +/- 2%: +/-0.02")
lo, hi = ratio(zs_of(z1, C1, 2.0)), ratio(zs_of(z1, C1, 0.5))
out += ["", f"combined: the payday convention alone spans ratio {min(lo,hi):.2f} - {max(lo,hi):.2f}; Claim 1 exponent (0.41-0.85) spans "
        f"{ratio(base, beta=0.41):.2f} - {ratio(base, beta=0.85):.2f}"]
solid = (max(lo, hi)/min(lo, hi) < 1.24)
out.append(f"DECISION: the 24% is {'SOLID' if solid else 'NOT solid'} -> {'run' if solid else 'do NOT run'} the today's-data test.")
txt = "\n".join(out); print(txt); open("iter63_audit.txt", "w").write(txt + "\n")
