"""Which dark-energy form can live in AeST's action? AeST's aether sees the expansion rate through Q = ∇·A = 3H.
 (a) ρ_DE ∝ ȧ^-β = (aH)^-β : our law (needs a as well: in AeST, a ∝ J^-1/3 with J the scalar condensate's conserved charge)
 (b) ρ_DE ∝ H^-β           : aether-only term F(Q)
Scored on DESI DR2 BAO + DES-Dovekie SN + CMB acoustic angle θ*, H0 fitted, ω_b and ω_c fixed from the CMB."""
import numpy as np
from scipy.optimize import minimize_scalar, brentq
import sorkin_test as S
def H_form(h, beta, form):
    a = np.exp(S.LNA); om, orad = (S.wb+S.wc)/h**2, S.wr/h**2; ol = 1-om-orad; E = np.empty_like(a)
    for i, ai in enumerate(a):
        m = om*ai**-3 + orad*ai**-4
        if form == "adot": f = lambda e: e*e - m - ol*(ai*e)**-beta
        else:              f = lambda e: e*e - m - ol*e**-beta
        E[i] = brentq(f, 1e-3, 1e20, rtol=1e-12)
    return h*100e3/S.Mpc*E
def best(beta, form):
    r = minimize_scalar(lambda h: S.score(S.LNA, H_form(h, beta, form))["tot"], bounds=(0.6, 0.75), method="bounded", options={"xatol": 1e-4})
    return S.score(S.LNA, H_form(r.x, beta, form))
ref = best(0.0, "H")["tot"]
print(f"ΛCDM χ² = {ref:.2f}")
for form, lab in [("adot", "ρ ∝ ȧ^-β (our law)"), ("H", "ρ ∝ H^-β (aether only)")]:
    rows = []
    for beta in [0.25, 0.5, 0.75, 1.0, 1.5]:
        s = best(beta, form); rows.append((beta, s["tot"]-ref, s["H0"]))
    bb = min(rows, key=lambda x: x[1])
    print(lab + ":  " + "  ".join(f"β={b}: Δχ² {d:+.1f} (H0 {h:.1f})" for b, d, h in rows))
