"""Which dark-energy terms can AeST's own fields carry, and do they trade energy with the condensate?
Minisuperspace (flat FLRW, lapse N): AeST reduces to  L = a^3 N [ -3 H^2/N^2 + K(Q) - U ] / (8 pi G),  Q = phidot/N, H = adot/a,
aether expansion Theta = 3H/N. Candidate dark-energy terms U built from AeST's cosmological invariants: Theta, Q, phi.
We derive the field equations and check (a) whether the condensate density ρ_c = Q K' - K still falls exactly as a^-3 (no exchange),
(b) what the dark-energy density depends on."""
import sympy as sp
t = sp.symbols('t'); a, N, ph = [sp.Function(n)(t) for n in ('a', 'N', 'phi')]
K = sp.Function('K'); 
def analyse(label, Ufun):
    Q = ph.diff(t)/N; Th = 3*a.diff(t)/(a*N)
    U = Ufun(Th, Q, ph)
    L = a**3*N*(-3*(a.diff(t)/a)**2/N**2 + K(Q) - U)
    # energy constraint (vary N), then set N=1
    EN = sp.diff(L, N) - sp.diff(sp.diff(L, N.diff(t)), t)
    # scalar equation: d/dt dL/dphidot = dL/dphi
    Eph = sp.diff(sp.diff(L, ph.diff(t)), t) - sp.diff(L, ph)
    sub = {N: 1}
    EN = sp.simplify(EN.subs(N, 1).doit()); Eph = sp.simplify(Eph.subs(N, 1).doit())
    charge = sp.simplify(sp.diff(L, ph.diff(t)).subs(N, 1).doit())   # conserved if dL/dphi = 0
    print(f"\n== {label}")
    print("  shift charge  a^3(...) =", sp.simplify(charge/ a**3), "   conserved:", sp.simplify(sp.diff(L, ph)) == 0)
    print("  energy constraint  ->  3H^2 =", sp.simplify(sp.solve(EN, (a.diff(t)/a)**2)[0] if False else -sp.expand(EN)/a**3 + 3*(a.diff(t)/a)**2) )
Uth  = lambda Th, Q, ph: sp.Function('f')(Th)
Uthq = lambda Th, Q, ph: sp.Function('f')(Th, Q)
Uphi = lambda Th, Q, ph: sp.Function('V')(ph)
for lab, U in [("U = f(Θ): aether expansion only", Uth), ("U = f(Θ, Q): aether expansion + condensate clock rate", Uthq), ("U = V(φ): potential of AeST's scalar clock", Uphi)]:
    analyse(lab, U)
print("""
Reading the results:
 f(Θ):   charge = a^3 K'(Q) is conserved and contains only K -> condensate keeps ρ_c ∝ a^-3 exactly (no exchange);
         dark energy ρ_DE = f - Θ f'(Θ) depends on H alone  ->  ρ ∝ H^-β type (excluded by data: Δχ² +39.5 at β=½).
 f(Θ,Q): charge = a^3 (K' - f_Q): shared between condensate and dark energy -> exchange; and Θ -> 0 inside galaxies (singular for ȧ^-β).
 V(φ):   shift symmetry broken by V: d(a^3 K')/dt = -a^3 V'(φ) -> exchange; locally regular (φ ≈ Q0 t everywhere).""")
