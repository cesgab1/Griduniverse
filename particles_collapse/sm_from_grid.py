"""
What can the grid fix about the Standard Model's free numbers?
 A. Charges: grid consistency (gauge-anomaly cancellation, needed for chiral fermions on a grid) + Higgs Yukawas -> hypercharges.
 B. Coupling strengths: run g_Y, g_2, g_3 (two-loop SM) up to the grid scale 1/l_d; test 'one coupling at the grid scale'.
 C. Higgs mass: grid boundary condition  lambda(mu_grid) = 0  (flat Higgs potential at the grid scale), and the stronger
    'multiple-point' version lambda = beta_lambda = 0. Two-loop running from NNLO MS-bar inputs (Buttazzo et al. 2013 fits).
 D. Generations: hard counting constraints.
Inputs (PDG 2024): m_t = 172.57 +- 0.29 GeV, m_H = 125.20 +- 0.11 GeV, alpha_s(M_Z) = 0.1180 +- 0.0009, M_W = 80.369 GeV.
"""
import numpy as np, sympy as sp, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

print("=== A. Charges from grid consistency ===")
yQ, yu, yd, yL, ye, yH = sp.symbols("yQ yu yd yL ye yH")
eqs = [2*yQ + yu + yd,                                   # [SU(3)]^2 U(1)
       3*yQ + yL,                                        # [SU(2)]^2 U(1)
       6*yQ + 3*yu + 3*yd + 2*yL + ye,                   # gravity^2 U(1)
       6*yQ**3 + 3*yu**3 + 3*yd**3 + 2*yL**3 + ye**3]    # [U(1)]^3
sols = sp.solve(eqs, [yu, yd, yL, ye], dict=True)
print("anomaly-free solutions (yQ free):", sols)
# Higgs Yukawas: Q H d^c, Q H~ u^c, L H e^c must be neutral
yuk = [yQ + yH + yd, yQ - yH + yu, yL + yH + ye]
for s in sols:
    sub = [sp.simplify(e.subs(s)) for e in yuk]
    solH = sp.solve(sub, [yH], dict=True)
    ok = [x for x in solH if all(sp.simplify(e.subs(x)) == 0 for e in sub)]
    print(f"  {s}: Yukawas allowed -> {'yes, Higgs yH = ' + str(ok[0][yH]) if ok else 'no'}")
s = [x for x in sols if x[yu] == -4*yQ][0]
scale = sp.Rational(1, 6) / yQ                            # normalise: electric charge Q = T3 + Y, electron charge -1
Yv = {k: sp.nsimplify(v.subs(yQ, sp.Rational(1, 6))) for k, v in s.items()}
qu = sp.Rational(1, 2) + sp.Rational(1, 6); qd = -sp.Rational(1, 2) + sp.Rational(1, 6); qe = -Yv[ye]; qnu = sp.Rational(1, 2) + Yv[yL]
print(f"  => up quark {qu}, down quark {qd}, electron {qe}, neutrino {qnu};  proton = {2*qu + qd}, neutron = {qu + 2*qd}, hydrogen atom = {2*qu + qd + qe}")

print("\n=== B & C. Running couplings to the grid scale ===")
mt, mH, a3, MW = 172.57, 125.20, 0.1180, 80.369
def inputs(mt=mt, mH=mH, a3=a3):
    lam = 0.12604 + 0.00206*(mH - 125.15) - 0.00004*(mt - 173.34)
    yt = 0.93690 + 0.00556*(mt - 173.34) - 0.00042*(a3 - 0.1184)/0.0007
    g3 = 1.1666 + 0.00314*(a3 - 0.1184)/0.0007 - 0.00046*(mt - 173.34)
    g2 = 0.64779 + 0.00004*(mt - 173.34) + 0.00011*(MW - 80.384)/0.014
    gY = 0.35830 + 0.00011*(mt - 173.34) - 0.00020*(MW - 80.384)/0.014
    return [gY, g2, g3, yt, lam]
k = 1/(16*np.pi**2)
def rge(t, y):
    gY, g2, g3, yt, l = y; Y2, G2, S2, T2 = gY*gY, g2*g2, g3*g3, yt*yt
    b1 = [41/6*gY**3, -19/6*g2**3, -7*g3**3,
          yt*(9/2*T2 - 8*S2 - 9/4*G2 - 17/12*Y2),
          24*l*l - 6*T2*T2 + 3/8*(2*G2*G2 + (G2 + Y2)**2) + l*(12*T2 - 9*G2 - 3*Y2)]
    b2 = [gY**3*(199/18*Y2 + 9/2*G2 + 44/3*S2 - 17/6*T2),
          g2**3*(3/2*Y2 + 35/6*G2 + 12*S2 - 3/2*T2),
          g3**3*(11/6*Y2 + 9/2*G2 - 26*S2 - 2*T2),
          yt*(-12*T2*T2 + T2*(131/16*Y2 + 225/16*G2 + 36*S2 - 12*l) + 1187/216*Y2*Y2 - 23/4*G2*G2 - 108*S2*S2
              - 3/4*Y2*G2 + 9*G2*S2 + 19/9*Y2*S2 + 6*l*l),
          (-312*l**3 - 144*l*l*T2 + 36*l*l*(3*G2 + Y2) - 3*l*T2*T2 + 80*l*S2*T2 + 45/2*l*G2*T2 + 85/6*l*Y2*T2
           - 73/8*l*G2*G2 + 39/4*l*G2*Y2 + 629/24*l*Y2*Y2 + 30*T2**3 - 32*S2*T2*T2 - 8/3*Y2*T2*T2 - 9/4*G2*G2*T2
           + 21/2*G2*Y2*T2 - 19/4*Y2*Y2*T2 + 305/16*G2**3 - 289/48*G2*G2*Y2 - 559/48*G2*Y2*Y2 - 379/48*Y2**3)]
    return [k*a + k*k*b for a, b in zip(b1, b2)]
EP = 1.220890e19; mu_grid = EP / 0.603
def run(y0, mu_to, mu_from=None, dense=False):
    mu_from = mu_from or 173.34
    return solve_ivp(rge, [np.log(mu_from), np.log(mu_to)], y0, rtol=1e-10, atol=1e-12, dense_output=dense)
sol = run(inputs(), mu_grid, dense=True)
tt = np.linspace(np.log(173.34), np.log(mu_grid), 3000); Y = sol.sol(tt)
lam = Y[4]; i0 = np.where(lam < 0)[0]
print(f"validation: lambda turns negative at {np.exp(tt[i0[0]]):.1e} GeV (literature ~1e10-1e11), lambda(M_Planck) = {np.interp(np.log(EP), tt, lam):+.4f} (literature ~ -0.013)")
gY, g2, g3 = Y[0, -1], Y[1, -1], Y[2, -1]
a1 = 5/3*gY**2/(4*np.pi); a2 = g2**2/(4*np.pi); a3g = g3**2/(4*np.pi)
aem = (gY*g2)**2/(gY**2 + g2**2)/(4*np.pi)
print(f"at the grid scale {mu_grid:.2e} GeV:  1/alpha_1(GUT-normalised) = {1/a1:.1f},  1/alpha_2 = {1/a2:.1f},  1/alpha_3 = {1/a3g:.1f},  1/alpha_em = {1/aem:.1f}")
spread = (max(1/a1, 1/a2, 1/a3g) - min(1/a1, 1/a2, 1/a3g)) / np.mean([1/a1, 1/a2, 1/a3g])
print(f"'one coupling at the grid scale': the three differ by {spread*100:.0f}% (would need 0%) -> fails as for plain SM grand unification")
# where do pairs meet?
for (i, j, name) in [(0, 1, "alpha_1 = alpha_2"), (1, 2, "alpha_2 = alpha_3"), (0, 2, "alpha_1 = alpha_3")]:
    ai = lambda t: (5/3 if i == 0 else 1)*sol.sol(t)[i]**2; aj = lambda t: (5/3 if j == 0 else 1)*sol.sol(t)[j]**2
    try:
        tc = brentq(lambda t: ai(t) - aj(t), np.log(200), np.log(mu_grid)); print(f"  {name} at {np.exp(tc):.1e} GeV")
    except ValueError: print(f"  {name}: never below the grid scale")

def lam_grid(mH_, mt_=mt, a3_=a3, mu=mu_grid):
    return run(inputs(mt_, mH_, a3_), mu).y[4, -1]
mH_pred = brentq(lambda m: lam_grid(m), 110, 150)
dm_t = (brentq(lambda m: lam_grid(m, mt + 0.29), 110, 150) - brentq(lambda m: lam_grid(m, mt - 0.29), 110, 150)) / 2
dm_a = (brentq(lambda m: lam_grid(m, mt, a3 + 0.0009), 110, 150) - brentq(lambda m: lam_grid(m, mt, a3 - 0.0009), 110, 150)) / 2
mt_pred = brentq(lambda m: lam_grid(125.20, m), 160, 185)
print(f"\nGrid boundary condition lambda(grid scale) = 0:")
print(f"  predicted Higgs mass = {mH_pred:.1f} ± {np.hypot(dm_t, dm_a):.1f} (inputs) ± ~1 (two-loop truncation) GeV   vs measured 125.20 ± 0.11")
print(f"  equivalently, predicted top-quark mass = {mt_pred:.2f} GeV   vs measured 172.57 ± 0.29 (and ~0.5-1 GeV pole-mass ambiguity)")
# multiple-point: lambda = 0 and beta_lambda = 0 at the same scale (scale free)
def mp(mH_):
    s_ = run(inputs(mt, mH_), mu_grid * 10, dense=True); t_ = np.linspace(np.log(173.34), np.log(mu_grid * 10), 4000); l_ = s_.sol(t_)[4]
    return l_.min(), np.exp(t_[np.argmin(l_)])
mH_mp = brentq(lambda m: mp(m)[0], 110, 150)
print(f"Stronger 'multiple-point' condition (lambda = beta_lambda = 0): m_H = {mH_mp:.1f} GeV, reached at {mp(mH_mp)[1]:.1e} GeV")

print("\n=== D. Number of generations: hard constraints ===")
for N in range(1, 9):
    cp = (N - 1)*(N - 2)//2; b0 = 11 - 4*N/3
    print(f"  N = {N}: CP-violating phases in quark mixing = {cp};  QCD asymptotically free: {'yes' if b0 > 0 else 'NO'} (b0 = {b0:.2f});  anomaly-free: yes (per generation)")
print("  measured: light neutrino species from the Z width = 2.984 ± 0.008 (LEP); a 4th chiral generation is excluded by the Higgs production rate")
json.dump(dict(charges=dict(u=str(qu), d=str(qd), e=str(qe)), inv_alpha_grid=[1/a1, 1/a2, 1/a3g, 1/aem], mH_pred=mH_pred, mH_err=float(np.hypot(dm_t, dm_a)),
               mt_pred=mt_pred, mH_mp=mH_mp), open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/sm_from_grid.json", "w"), indent=1)
