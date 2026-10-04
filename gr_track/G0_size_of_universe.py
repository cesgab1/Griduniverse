"""GR track G0 (Coalesce): how big does GR, with today's data, say the universe is? GR only -- no grid input."""
import numpy as np
from scipy.integrate import quad
c = 299792.458; H0 = 67.4; Om, Or = 0.315, 9.1e-5; OL = 1 - Om - Or
E = lambda a: np.sqrt(Om/a**3 + Or/a**4 + OL)
Mpc_Gly = 3.2616e-3
chi_p = c/H0*quad(lambda a: 1/(a*a*E(a)), 1e-12, 1, limit=400)[0]                 # particle horizon (comoving, Mpc)
chi_ls = c/H0*quad(lambda a: 1/(a*a*E(a)), 1/1090, 1, limit=400)[0]               # last scattering
eh = c/H0*quad(lambda x: 1/(np.exp(x)*E(np.exp(x))), 0, 60, limit=400)[0]                       # event horizon (comoving today)
age = quad(lambda a: 1/(a*E(a)), 1e-12, 1, limit=400)[0]/(H0/977.8)                # Gyr
out = ["GR track G0: the size of the universe according to GR + today's data", "",
       f"age: {age:.2f} billion years",
       f"observable universe (particle horizon): radius {chi_p/1000:.1f} Gpc = {chi_p*Mpc_Gly:.1f} billion light years (diameter {2*chi_p*Mpc_Gly:.0f})",
       f"oldest light (last scattering): radius {chi_ls/1000:.2f} Gpc",
       f"region we can EVER reach from now (event horizon): {eh/1000:.1f} Gpc = {eh*Mpc_Gly:.1f} billion light years",
       f"visible volume: {4/3*np.pi*(chi_p/1000)**3:.0f} Gpc^3", ""]
# curvature: DESI DR2 + CMB Omega_k = 0.0023 +/- 0.0011 (positive = open/negative curvature in this convention)
ok, sk = 0.0023, 0.0011
out.append("whole universe (GR's equations allow three shapes; data decide between them via curvature Omega_k = 0.0023 +/- 0.0011):")
for nsig in (2, 3):
    lo = ok - nsig*sk
    if lo < 0:
        R = c/H0/np.sqrt(-lo)/1000; V = 2*np.pi**2*R**3
        out.append(f"  {nsig} sigma: a closed (finite) universe is still allowed if its curvature radius > {R:.0f} Gpc -> total volume > {V:.2e} Gpc^3"
                   f" = {V/(4/3*np.pi*(chi_p/1000)**3):.0f} x the visible volume")
    else:
        out.append(f"  {nsig} sigma: closed shape EXCLUDED -> flat or open -> infinite, unless space wraps around (topology: GR is silent)")
out += ["", "where GR stops: (1) at time zero (the singularity -- infinite density, GR breaks down); (2) GR has no rule for the",
        "universe's global shape (topology) -- whether flat space goes on forever or wraps around is outside GR."]
txt = "\n".join(out); print(txt); open("G0_size_of_universe.txt", "w").write(txt + "\n")
