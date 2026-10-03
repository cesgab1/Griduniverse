"""
ITERATION 4: (A) what sets the memory?  (B) does the jostling relic need a TRUE thermal spectrum?  (C) mass and fifth force.

A. If the jostled tension is the VELOCITY of a grid field (tau = dPhi/dt), the expansion damps it by Hubble friction:
   minimally coupled: d/dt(a^3 Phi_dot) = force  ->  tau ∝ a^-3 without force: rate 3H  (kappa = 3)
   conformally coupled homogeneous mode: guessed rate 2H; the numerical check below gives ~0.9H instead (guess was wrong)
   Checked numerically below by integrating the field equations in a matter + Lambda background.
B. Iteration 2's random walk came from the k -> 0 thermal tail n(k) ≈ Theta/k (Rayleigh-Jeans). A relic made by decays
   without self-interaction has n(k) -> const at small k: test whether it still gives Var ∝ memory.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
out = ["ITERATION 4", "", "A. Hubble friction of a field velocity (no force), matter + Lambda background, from a = 0.1 to 1:"]
Om = 0.31; H = lambda a: np.sqrt(Om*a**-3 + 1 - Om)
for lab, xi in (("minimal coupling", 0.0), ("conformal coupling (xi = 1/6)", 1/6)):
    # homogeneous mode: Phi'' + 3H Phi' + xi R Phi = 0, R = 6(2H^2 + Hdot) ; use ln a as time variable
    def rhs(lna, y):
        a = np.exp(lna); h = H(a); dh = -1.5*Om*a**-3/h                         # dH/dt / H  -> Hdot = h*dh*? (dH/dlna = dh/h*h?)
        Hdot = -1.5*Om*a**-3                                                     # Hdot = dH/dt in units H0^2 for matter + Lambda
        R = 6*(2*h*h + Hdot); P, V = y                                           # V = dPhi/dt
        return [V/h, (-3*h*V - xi*R*P)/h]
    s = solve_ivp(rhs, [np.log(0.1), 0], [1.0, 1.0], rtol=1e-10, atol=1e-12, dense_output=True)
    la = np.linspace(np.log(0.3), 0, 50); V = np.abs(s.sol(la)[1]) + 1e-30
    rate = -np.gradient(np.log(V), la)                                           # = (decay rate of the velocity)/H
    out.append(f"   {lab:32s}: decay rate of Phi_dot / H, late times: {np.median(rate):.2f}")
out.append("   -> the memory is NOT a free choice once the tension is a field velocity: kappa = 3 (minimal coupling); a conformally")
out.append("      coupled tension field would give ~0.9, which the data disfavour (iteration 1: memory ~1 Hubble time barely beats Lambda).")
out.append("")
out.append("   Fits with kappa = 3 (iter1_fits_*_k3.json), Delta chi2 vs Lambda, NO free dark-energy parameter:")
import json
for sn in ("PANTHEON", "DESY5", "UNION3"):
    d = json.load(open(f"iter1_fits_{sn}_k3.json")); L = d["LCDM"]
    out.append(f"     {sn:9s}: linear energy {d['lin kappa=3.0'][0] - L:+.1f}   quadratic energy {d['quad kappa=3.0'][0] - L:+.1f}")
out.append("     (instant Claim 1: -5.4 / -7.2 / -6.8; free w0wa with 2 parameters: -7.1 / -10.3 / -13.2)")
out.append("")
out.append("B. Does the relic need a thermal (Rayleigh-Jeans) low-k tail?  Var vs memory, Theta = 1, R = 1:")
def var(gamma, occ):
    f = lambda k: (k/(4*np.pi**2))*2*occ(k)*np.exp(-k*k)/(gamma**2 + k*k)
    return quad(f, 0, 40, limit=500, points=[gamma, 1.0])[0]
thermal = lambda k: 1/np.expm1(k)                       # Bose-Einstein, ≈ 1/k at small k
decay = lambda k: 0.5*np.exp(-k)                         # free-streaming, made by decays: finite occupation at small k
for g in (1e-2, 1e-3, 1e-4, 1e-5):
    out.append(f"   memory {1/g:6.0e}:  thermal relic Var = {var(g, thermal):10.3e}    decay-made relic (no self-interaction) Var = {var(g, decay):.3e}")
out.append("   -> only a THERMAL relic (or one with a Rayleigh-Jeans tail, i.e. self-interacting) gives the random walk.")
out.append("")
out.append("C. Consequences for the relic scalar s:")
out.append("   * Inflation (required by the model) dilutes any earlier relic to nothing, and a conformal field is not made by")
out.append("     inflation itself. So s must be made at reheating (inflaton -> s) and THERMALISE itself (e.g. a lambda_s s^4 self-")
out.append("     coupling, which is conformal). Its temperature is then set by the inflaton's branching ratio, NOT by contact with")
out.append("     ordinary matter: iteration 3's 'Delta N_eff = 0.027-0.057' holds only if s also shared a temperature with ordinary")
out.append(f"     matter. In general 0 < Delta N_eff < 0.107, i.e. T_s/T_nu < {(0.107*7/4)**0.25:.2f}. The amplitude change is absorbed in")
out.append("     the (already free) size of dark energy. PREDICTION WEAKENED: Delta N_eff > 0, size not fixed.")
out.append("   * Fifth force: with a symmetry s -> -s (and conformal coupling to curvature), s has no LINEAR coupling to matter, so")
out.append("     no long-range fifth force at tree level. Contact with ordinary matter via a Higgs portal lambda s^2|H|^2 would give s")
out.append("     a mass ~ sqrt(lambda) x 246 GeV; staying massless enough for horizon-scale noise (m_s << kappa H0 ~ 1e-32 eV) needs")
out.append(f"     lambda << {(1e-32/246e9)**2:.0e}, far too weak ever to thermalise. So: NO Higgs-portal contact -> s is a separate")
out.append("     dark radiation, consistent with the general 0 < Delta N_eff < 0.107.")
out.append("   * Naturalness: a scalar lighter than 1e-32 eV is as unprotected as a quintessence field (same open problem).")
txt = "\n".join(out); print(txt); open("iter4_memory_and_source.txt", "w").write(txt + "\n")
