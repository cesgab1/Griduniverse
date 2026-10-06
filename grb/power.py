"""Can gamma-ray-burst distances tell our dark-energy law from constant dark energy? (power check before hunting for data)
Both models fitted to the same BAO+SN+CMB data (best fits from fit_law.py); compare distance-modulus SHAPE at z = 1-9
(an overall offset is absorbed by burst calibration)."""
import os, numpy as np
from scipy.integrate import cumulative_trapezoid
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../cosmology_fits/fit_law.py")).read()
exec(src.split("def chi2")[0])
def mu(model, Om, h):
    E = E_of_z(model, Om, h, []); chi = cumulative_trapezoid(1/E, ZG, initial=0)
    return lambda z: 5*np.log10((1 + z)*np.interp(z, ZG, chi))
mL, mW = mu("LCDM", 0.302, 0.684), mu("LAW", 0.310, 0.675)      # iteration 92 combined best fits (Pantheon+)
zs = np.array([1, 2, 3, 4, 6, 8, 9.0])
d = mW(zs) - mL(zs); d -= d[0]                                    # remove offset (calibration), anchor at z = 1
for z, x in zip(zs, d): print(f"z = {z:3.0f}: law minus constant dark energy, distance modulus (shape) {x:+.4f} mag")
print("burst distances: scatter ~0.5-1.0 mag per burst after calibration (energy-peak relations, mem);")
for N in (118, 220, 1000):
    print(f"  {N} bursts at z > 1: mean shape error ~ {0.8/np.sqrt(N/2):.3f} mag per half-sample")
