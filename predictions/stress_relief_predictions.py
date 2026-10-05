"""Predictions of the global stress-relief law d ln rho_DE / d ln a = beta q (Claim 1). Registered Oct 4 2026."""
import numpy as np
from scipy.integrate import solve_ivp
Om, h = 0.311, 0.674; Or = 9.1e-5; ODE = 1 - Om - Or
def hist(beta):
    sol = solve_ivp(lambda x, y: [beta*((Om*np.exp(-3*x)/2 + Or*np.exp(-4*x) - np.exp(y[0]))/(Om*np.exp(-3*x) + Or*np.exp(-4*x) + np.exp(y[0])))],
                    [0, np.log(1/6)], [np.log(ODE)], dense_output=True, rtol=1e-10, atol=1e-12)
    z = np.linspace(0, 5, 5001); x = -np.log(1+z)
    rde = np.exp(sol.sol(x)[0]); rm = Om*(1+z)**3; rr = Or*(1+z)**4
    E = np.sqrt(rm + rr + rde)
    q = (rm/2 + rr - rde*(1 + 1.5*(-beta*((rm/2 + rr - rde)/(rm+rr+rde))/3)))/(E**2)   # q incl. DE pressure
    w = -1 - beta*((rm/2 + rr - rde)/(rm + rr + rde))/3
    return z, E, rde/ODE, w, q
out = ["Predictions of the stress-relief (Claim 1) law, registered 2026-10-04 before DESI DR3 / Euclid / ELT results.", ""]
for beta in (0.35, 0.5, 0.85):
    z, E, r, w, q = hist(beta)
    iw = np.argmin(abs(w + 1)); iq = np.argmin(abs(q))
    out.append(f"beta = {beta}: w today {w[0]:.3f}; w = -1 crossing at z = {z[iw]:.3f}; acceleration starts (q = 0) at z = {z[iq]:.3f}; "
               f"dark-energy peak {100*(r.max()-1):.1f}% above today at z = {z[np.argmax(r)]:.2f}; w(z=1) = {np.interp(1, z, w):.3f}")
# redshift drift over 20 years, Claim 1 (beta 0.5) vs constant
z, E, *_ = hist(0.5); zL = z; EL = np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + ODE)
H0s = h*100/3.0857e19; yr = 3.156e7; c = 2.998e8
for zz in (1, 2, 3, 4):
    dv = lambda EE: c*H0s*20*yr*(1 - np.interp(zz, z, EE)/(1+zz))
    out.append(f"redshift drift over 20 yr at z = {zz}: Claim 1 {100*dv(E):.2f} cm/s, constant {100*dv(EL):.2f} cm/s")
out.append("")
out.append("Fixed features (any beta > 0): (1) dark energy's w crosses -1 at EXACTLY the redshift where acceleration begins (w = -1 - beta q/3);")
out.append("(1b) dark energy PEAKS at that same redshift; (2) no Big Rip, no crunch (beta = 1/2: far future a ~ t^5); (3) open curvature Omega_k ~ +0.002 (P1);")
out.append("(4) gravitational waves at c with only GR's two polarisations (passed); (5) dark energy uniform -- no extra push near masses.")
txt = "\n".join(out); print(txt); open("stress_relief_predictions.txt", "w").write(txt + "\n")
