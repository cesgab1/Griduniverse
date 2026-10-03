"""
Can the Tension-Rate Law come from an ACTION with a preferred time slicing (Hořava / khronometric / VCDM-type)?
Minisuperspace (flat FRW): lapse N(t), scale factor a(t), expansion of the time slices K = 3H, H = adot/(N a).
Grid ingredient: physical link length ell = ell0 * a^s (s = 1: links stretch with space; s = 0: cells are added).
  S = ∫ dt N a^3 [ -3 M^2 H^2 + P(H, a) ] - ∫ dt N a^3 rho_m(a)
Checks (sympy):
 1. Vary N -> Friedmann constraint 3 M^2 H^2 = rho_m + rho_DE with rho_DE = H dP/dH - P.
 2. P = -(2/3) C (H a^s)^(-1/2) gives rho_DE = C (H a^s)^(-1/2): for s = 1, rho ∝ adot^(-1/2) (Claim 1). Unique up to a
    total derivative.
 3. Vary a -> second equation; check it is CONSISTENT with the constraint (constraint preserved in time) and gives
    w_DE = -1 - (q + 1 - s)/6 (Claim 1 for s = 1).
 4. Analytic EFT: P = sum c_n H^n (n >= 0) gives rho_DE = sum (n-1) c_n H^n: only H^0, H^2, H^3...; never H^(-1/2).
"""
import sympy as sp
t = sp.symbols("t"); M, C, s, rho0 = sp.symbols("M C s rho0", positive=True)
N = sp.Function("N")(t); a = sp.Function("a")(t)
Hs, As = sp.symbols("H A", positive=True)
P = -sp.Rational(2, 3)*C*(Hs*As**s)**sp.Rational(-1, 2)
rhoDE_formula = sp.simplify(Hs*sp.diff(P, Hs) - P)
print("1-2. rho_DE = H dP/dH - P =", rhoDE_formula, "   (target C (H a^s)^(-1/2))")
H = a.diff(t)/(N*a)
rho_m = rho0*a**-3
Lag = N*a**3*(-3*M**2*H**2 + P.subs({Hs: H, As: a})) - N*a**3*rho_m
EN = sp.diff(Lag, N) - sp.diff(sp.diff(Lag, N.diff(t)), t)
Ea = sp.diff(Lag, a) - sp.diff(sp.diff(Lag, a.diff(t)), t)
gauge = {N: 1}
constraint = sp.simplify(EN.subs(N, 1).doit())
# constraint should be: -(3 M^2 H^2 - rho_m - rho_DE) * a^3
h = sp.symbols("h", positive=True)
c_check = sp.simplify(sp.powsimp(sp.expand((constraint/a**3).subs(a.diff(t), h*a) - (3*M**2*h**2 - rho_m - C*(h*a**s)**sp.Rational(-1, 2))).subs(a, As), force=True))
print("1. constraint matches Friedmann with rho_DE above:", c_check == 0)
# 3. Noether identity: d/dt(constraint) must be a combination of the a-equation -> constraint preserved
Ea1 = sp.simplify(Ea.subs(N, 1).doit())
ident = sp.simplify(sp.powsimp(sp.expand(sp.diff(constraint, t) - a.diff(t)*Ea1), force=True))      # reparametrisation identity: dE_N/dt = adot * E_a (N=1)
print("3. Bianchi/Noether identity  d(constraint)/dt - adot * (a-equation) =", ident)
# effective pressure: a-equation = -3 a^2 (2 M^2 (addot/a) + M^2 H^2 ... ) form; extract w_DE from continuity instead:
q = sp.symbols("q"); lnrho = sp.Rational(-1, 2)*(sp.log(Hs) + s*sp.log(As))
# d ln H / d ln a = -(1 + q)
dlnrho_dlna = sp.Rational(-1, 2)*(-(1 + q) + s)
w = sp.simplify(-1 - dlnrho_dlna/3)
print("3. w_DE =", w, "  -> s = 1:", sp.simplify(w.subs(s, 1)), " (Claim 1: -1 - q/6)")
# 4. analytic EFT
n = sp.symbols("n", integer=True, nonnegative=True); cn = sp.symbols("c_n")
print("4. P = c_n H^n  ->  rho =", sp.simplify(Hs*sp.diff(cn*Hs**n, Hs) - cn*Hs**n), " : powers n = 0, 2, 3, ... only (n = 1 is a total derivative)")
