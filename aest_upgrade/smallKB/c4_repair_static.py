"""
Repair attempt: add c4 J^mu J_mu to AeST (Einstein-aether sign: +c4|J|^2 in L). alpha1 = -4[K_B + c4 + (2-K_B)/(1+lambda_s)].
Static weak field (A^mu static => J_i = d_i Phi, Psi = Phi). Quasistatic Lagrangian x(-16 pi G):
   (2-K_B-c4)|grad Phi|^2 - 2(2-K_B) grad Phi.grad phi + (2-K_B)|grad phi|^2 + (2-K_B) Jfun(Y),  plus 16 pi G Phi rho
Field equations (spherical): grad phi (1+J') = grad Phi ;  div[ ((2-K_B-c4) - (2-K_B)/(1+J')) grad Phi ] = 8 pi G rho
c4 = 0 reproduces AeST: coefficient (2-K_B) J'/(1+J') -> deep MOND (J' ~ |grad phi|/a0).
"""
import sympy as sp
KB, c4, Jp = sp.symbols("K_B c_4 Jp", real=True)
coef = (2-KB-c4) - (2-KB)/(1+Jp)
print("coefficient of grad Phi:", sp.simplify(coef))
print("deep-MOND limit (J' -> 0):", sp.limit(coef, Jp, 0), " -> must vanish for MOND; with c4 = -K_B it is", sp.limit(coef, Jp, 0).subs(c4, -KB))
print("\nConsequence: residual Newtonian term -c4 makes gravity Newtonian with G' = 2G/|c4| below g ~ |c4| a0/(2-K_B).")
for kb in (0.5, 1e-3, 2.5e-5):
    print(f"   c4 = -K_B = {-kb:g}: MOND lost below g ~ {kb/(2-kb):.1e} a0, replaced by Newton x {2/kb:.0f}")
print("\nalpha1 = -4(K_B + c4) (screened). Keeping MOND needs |c4| << 1 -> alpha1 ~ -4 K_B again -> K_B < 2.5e-5 -> cosmological")
print("instability (growth_exponent.txt). The c4 repair does not escape: MOND, alpha1 and early-universe stability conflict.")
