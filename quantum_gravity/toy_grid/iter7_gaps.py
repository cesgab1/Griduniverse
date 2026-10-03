"""
ITERATION 7 (Coalesce): the gap between layers is big enough that each layer barely interferes with the ones above or below.
Iteration 6 required correlation c <~ 6e-6 assuming EVERY pair of layers shares the same correlation (a common jostle).
With gaps, correlation instead FALLS OFF with distance in the stack: c_k between layers k apart (e.g. c_k = r^k, leakage
fraction r per gap). General result (Gaussian tension, quadratic energy):
    delta rho/rho = sqrt( 2 [1 + 2 sum_{k>=1} c_k^2 (1 - k/N)] / N )
-> only the total neighbour leakage matters; the effective number of independent layers is N / (1 + 2 sum c_k^2).
Monte Carlo: AR(1) chain across the stack (each layer = r x the layer below + fresh jostle).
"""
import numpy as np
rng = np.random.default_rng(5); out = ["ITERATION 7: gaps between layers (Coalesce)", ""]
def mc(N, r, patches=40000):
    t = np.empty((patches, N)); t[:, 0] = rng.standard_normal(patches)
    for i in range(1, N): t[:, i] = r*t[:, i - 1] + np.sqrt(1 - r*r)*rng.standard_normal(patches)
    e = (t**2).sum(1); return e.std()/e.mean()
def formula(N, r):
    k = np.arange(1, N); return np.sqrt(2*(1 + 2*np.sum(r**(2*k)*(1 - k/N)))/N)
out.append(" layers N  leakage per gap r   Monte Carlo   formula   effective independent layers")
for N, r in ((1000, 0.0), (1000, 0.5), (1000, 0.9), (1000, 0.99), (2000, 0.9)):
    m, f = mc(N, r), formula(N, r); out.append(f"   {N:5d}        {r:4.2f}           {m:.4f}      {f:.4f}     {2/f**2:8.0f}")
out.append("")
out.append("Compare a COMMON jostle shared by all layers (iteration 6): correlation c for every pair -> lumps sqrt(2(1/N + c^2)):")
for c in (0.01, 0.1):
    out.append(f"   N = 1e10, common correlation {c}: lumps {np.sqrt(2*(1e-10 + c*c)):.2e}  (does not average away)")
for r in (0.5, 0.9, 0.99):
    ell = 1 + 2*r*r/(1 - r*r)
    out.append(f"   N = 1e10, neighbour leakage r = {r}: lumps {np.sqrt(2*ell/1e10):.2e}  (averages away; costs a factor {np.sqrt(ell):.1f})")
need = lambda r: 2.6e10*(1 + 2*r*r/(1 - r*r))
out.append("")
out.append("Layers needed (CMB imprint < 1e-6) vs leakage per gap: " + ", ".join(f"r={r}: {need(r):.1e}" for r in (0, 0.5, 0.9, 0.99)))
out.append("-> Gaps rescue the picture: neighbour leakage only multiplies the needed number of layers by (1 + 2 sum c_k^2).")
out.append("   What must be ABSENT is a jostle COMMON to all layers (one shared relic, or a shared long-range mode).")
txt = "\n".join(out); print(txt); open("iter7_gaps.txt", "w").write(txt + "\n")
