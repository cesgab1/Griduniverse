"""Upgraded AeST dark energy, written as a potential V(φ) of AeST's scalar clock (φ ≈ Q0 t everywhere, so no local singularity).
V is chosen so that on the cosmic solution ρ_DE ∝ ȧ^(-1/2). Because V(φ) has w = -1 exactly, energy conservation forces an exchange with the
dark-matter-like condensate:   d(ρ_c a^3)/da = -a^3 dρ_DE/da.
Compared on DESI DR2 BAO + DES-Dovekie SN + θ*: ΛCDM, the law as a separate fluid (our earlier fits), and the law with exchange (AeST form)."""
import numpy as np
from scipy.optimize import minimize_scalar, brentq
import sorkin_test as S
def H_exchange(h, beta):
    a = np.exp(S.LNA); ob, oc, orad = S.wb/h**2, S.wc/h**2, S.wr/h**2
    # iterate: guess ρ_DE(a) from the law with E(a), then ρ_c from exchange, then E(a); normalise ρ_DE today by flatness
    E = np.sqrt((ob+oc)*a**-3 + orad*a**-4 + (1-ob-oc-orad))
    for _ in range(60):
        adot = a*E; shape = adot**-beta; shape /= shape[-1]
        # ρ_c a^3 = oc(early) - ∫ a^3 dρ_DE ; ρ_DE = ode*shape, ode fixed by E(1)=1
        def build(ode):
            rde = ode*shape; da = np.diff(rde)
            integ = np.concatenate([[0], np.cumsum(0.5*(a[1:]**3 + a[:-1]**3)*da)])
            rc = (oc - integ)/a**3        # oc = early-universe value (CMB)
            return rde, rc
        ode = brentq(lambda o: ob + build(o)[1][-1] + orad + o - 1, 0.3, 1.0)
        rde, rc = build(ode)
        En = np.sqrt(ob*a**-3 + rc + orad*a**-4 + rde)
        if np.max(np.abs(En/E - 1)) < 1e-10: break
        E = En
    return h*100e3/S.Mpc*En, rc[-1]*h**2, ode
def H_fluid(h, beta):
    a = np.exp(S.LNA); om, orad = (S.wb+S.wc)/h**2, S.wr/h**2; ol = 1-om-orad
    E = np.array([brentq(lambda e: e*e - om*ai**-3 - orad*ai**-4 - ol*(ai*e)**-beta, 1e-3, 1e20, rtol=1e-12) for ai in a])
    return h*100e3/S.Mpc*E
def best(fn):
    r = minimize_scalar(lambda h: S.score(S.LNA, fn(h))["tot"], bounds=(0.6, 0.75), method="bounded", options={"xatol": 1e-4})
    return r.x, S.score(S.LNA, fn(r.x))
hL, L = best(lambda h: H_fluid(h, 0.0)); print(f"ΛCDM: χ² {L['tot']:.2f}  (BAO {L['bao']:.1f}, SN {L['sn']:.1f}) H0 {L['H0']:.2f}")
for beta in [0.5]:
    hF, F = best(lambda h: H_fluid(h, beta))
    hX, X = best(lambda h: H_exchange(h, beta)[0])
    _, wc0, ode = H_exchange(hX, beta)
    print(f"β={beta} separate fluid : Δχ² {F['tot']-L['tot']:+.1f} (BAO {F['bao']-L['bao']:+.1f}, SN {F['sn']-L['sn']:+.1f}) H0 {F['H0']:.2f}")
    print(f"β={beta} AeST exchange  : Δχ² {X['tot']-L['tot']:+.1f} (BAO {X['bao']-L['bao']:+.1f}, SN {X['sn']-L['sn']:+.1f}) H0 {X['H0']:.2f};"
          f"  condensate today ω_c = {wc0:.4f} vs early {S.wc:.4f} ({100*(wc0/S.wc-1):+.1f}%)")
np.save("aest_de_best.npy", np.array([hX]))
