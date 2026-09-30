"""Clumping in upgraded AeST: sub-horizon growth of condensate + baryons with the energy exchange (transfer along the condensate's
motion, homogeneous, as for a clock potential V(φ)):  δc' = -θc/H - (Q/ρc H) δc ;  no momentum transfer.  Same early amplitude as ΛCDM."""
import numpy as np, json
from scipy.integrate import solve_ivp
import sorkin_test as S, aest_de as A
hX = float(np.load("aest_de_best.npy")[0]); hL = 0.6807
def run(h, beta):
    a = np.exp(S.LNA); ob, orad = S.wb/h**2, S.wr/h**2
    if beta == 0:
        oc = S.wc/h**2; E = np.sqrt((ob+oc)*a**-3 + orad*a**-4 + 1-ob-oc-orad); rc = oc*a**-3; rde = (1-ob-oc-orad)*np.ones_like(a)
    else:
        H, _, ode = A.H_exchange(h, beta); E = H/(h*100e3/S.Mpc)
        adot = a*E; shape = adot**-beta; shape /= shape[-1]; rde = ode*shape
        rc = E**2 - ob*a**-3 - orad*a**-4 - rde
    lna = S.LNA; dlnrde = np.gradient(np.log(rde), lna)
    Qfac = -rde*dlnrde/rc            # Q/(ρc H): energy into condensate per Hubble time, in units of ρc
    Ef = lambda N: np.interp(N, lna, E); Oc = lambda N: np.interp(N, lna, rc)/Ef(N)**2; Ob = lambda N: ob*np.exp(-3*N)/Ef(N)**2
    q = lambda N: np.interp(N, lna, Qfac); dlnE = lambda N: np.interp(N, lna, np.gradient(np.log(E), lna))
    def f(N, y):
        dc, tc, db, tb = y     # t = θ/H
        src = 1.5*(Oc(N)*dc + Ob(N)*db)
        return [-tc - q(N)*dc, -tc*(2+dlnE(N)) - src, -tb, -tb*(2+dlnE(N)) - src]
    N0 = np.log(1/300.)
    s = solve_ivp(f, [N0, 0], [1/300., -1/300., 1/300., -1/300.], rtol=1e-9, atol=1e-14)
    dc, db = s.y[0, -1], s.y[2, -1]; Om0 = rc[-1] + ob
    return (rc[-1]*dc + ob*db)/Om0, Om0
DL, OmL = run(hL, 0.0); DX, OmX = run(hX, 0.5)
s8L = 0.811          # ΛCDM σ8 for this early amplitude (CAMB, earlier fits)
s8X = s8L*DX/DL
print(f"ΛCDM: σ8 {s8L:.3f}, Ω_m {OmL:.3f}, S8 {s8L*np.sqrt(OmL/0.3):.3f}")
print(f"upgraded AeST (β=½, exchange): growth ×{DX/DL:.3f} -> σ8 {s8X:.3f}, Ω_m {OmX:.3f}, S8 {s8X*np.sqrt(OmX/0.3):.3f}")
print("measured S8: KiDS-Legacy 0.815 (+0.016 -0.021); DES Y3 0.776 ± 0.017; eROSITA clusters 0.86 ± 0.01")
