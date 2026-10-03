"""
ITERATION 6 (Coalesce's feedback): each LAYER has its own jostle; stacked ('pancaked') together they give a smoother,
more harmonious wave. Escape (a) [many independent components] becomes escape (b) [a smooth global response].

A. N stacked layers, each with its own jostled tension tau_i (Gaussian, same variance), quadratic energy
   rho = sum_i tau_i^2. If the layers' jostles have pairwise correlation c, the patch-to-patch lumpiness is
       delta rho / rho = sqrt(2 (1/N + c^2 (N-1)/N))    (Gaussian: Cov(x^2, y^2) = 2 c^2; Monte Carlo check below)
   Single layer: sqrt(2) -> CMB imprint ~0.16 (iteration 5). Requirement: imprint < 1e-6 (10% of the observed 1e-5).
B. Each layer needs its own jostlers (its own relic radiation), otherwise their jostles are correlated (c = 1).
   N relic species: Delta N_eff = N (4/7) (T_s/T_nu)^4 < 0.107.
C. Escape (c): can short-wavelength jostling keep the 1/adot law AND have a short correlation length?
   Noise weight per wavenumber ∝ k^p (p = 0 is thermal Rayleigh-Jeans): measure how Var scales with the memory rate gamma
   (the law needs Var ∝ 1/gamma) and the correlation length.
"""
import numpy as np
from scipy.integrate import quad
rng = np.random.default_rng(3); out = ["ITERATION 6: stacked layers (Coalesce) and short-wave jostling", ""]
lump = lambda N, c: np.sqrt(2*(1/N + c*c*(N - 1)/N))
out.append("A. Monte Carlo check of the lumpiness formula (200,000 patches):")
for N, c in ((1, 0), (10, 0), (100, 0), (100, 0.01), (1000, 0.001)):
    z0 = rng.standard_normal(200000); tau = np.sqrt(c)*z0[:, None] + np.sqrt(1 - c)*rng.standard_normal((200000, N)) if c > 0 else rng.standard_normal((200000, N))
    r = (tau**2).sum(1); out.append(f"   N = {N:5d}, correlation {c:6.3f}: Monte Carlo {r.std()/r.mean():.4f}   formula {lump(N, c):.4f}")
Phi1 = 0.16
out.append("")
out.append(f"   CMB imprint = {Phi1} x lumpiness/sqrt(2). Requirement imprint < 1e-6:")
for N in (74, 1e4, 1e8, 1e10, 1e12):
    out.append(f"   N = {N:8.0e} independent layers: imprint {Phi1*lump(N, 0)/np.sqrt(2):.1e}   {'OK' if Phi1*lump(N, 0)/np.sqrt(2) < 1e-6 else 'too lumpy'}")
Nreq = (Phi1/1e-6)**2; out.append(f"   -> need N >~ {Nreq:.1e} layers, AND layer-to-layer correlation c <~ {1/np.sqrt(Nreq):.0e} (else the shared part stays lumpy)")
out.append("   (74 = the layer count suggested earlier by the electron-mass hierarchy: far too few, if those are the same layers)")
out.append("")
out.append("B. Each layer with its own relic radiation: Delta N_eff = N (4/7)(T_s/T_nu)^4 < 0.107 ->")
for N in (1e10, 1e12):
    out.append(f"   N = {N:.0e}: T_s/T_nu < {(0.107*7/4/N)**0.25:.1e}  (T_s < {(0.107*7/4/N)**0.25*1.95*1e3:.1f} mK): allowed; amplitude absorbed in the free height")
out.append("")
out.append("C. Short-wave jostling: noise weight ∝ k^p, memory filter 1/(gamma^2 + k^2), cell cutoff exp(-k^2):")
def var(g, p): return quad(lambda k: k**p*np.exp(-k*k)/(g*g + k*k), 0, 30, limit=600, points=[g])[0]
def corr_len(g, p):
    c0 = var(g, p); f = lambda d: quad(lambda k: k**p*np.exp(-k*k)*np.sin(k*d)/(k*d)/(g*g + k*k), 0, 30, limit=2000, points=[g])[0]/c0 - 0.5
    lo, hi = 1e-3, 1e4
    for _ in range(60):
        mid = np.sqrt(lo*hi); (lo, hi) = (mid, hi) if f(mid) > 0 else (lo, mid)
    return np.sqrt(lo*hi)
for p in (0.0, 0.5, 1.0, 1.5, 2.0):
    g1, g2 = 1e-3, 1e-4; ex = np.log(var(g2, p)/var(g1, p))/np.log(g1/g2)
    out.append(f"   p = {p:3.1f}: Var ∝ gamma^-{ex:.2f} (law needs 1)   half-correlation length at gamma=1e-3: {corr_len(g1, p):8.1f}  (c/gamma = 1000, cell = 1)")
out.append("   -> keeping the 1/adot law (exponent 1) forces a correlation length ~ c/gamma (horizon/3); short correlation")
out.append("      lengths lose the law (exponent -> 0: Lambda-like). Escape (c) is CLOSED.")
txt = "\n".join(out); print(txt); open("iter6_layers.txt", "w").write(txt + "\n")
