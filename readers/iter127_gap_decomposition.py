"""PREREG 127: decompose the gap between the measured CP bias and the needed matter excess."""
import numpy as np
from scipy.optimize import brentq
J = 3.18e-5
pole = dict(u=0.00216, d=0.00470, s=0.0935, c=1.2730, b=4.183, t=172.57)      # PDG 2024 (sourced)
run  = dict(u=0.00127, d=0.00271, s=0.0553, c=0.619, b=2.86, t=171.7)         # ~at 91 GeV, FROM MEMORY (approx.)
def D(T, m):
    u,d,s,c,b,t = (m[k] for k in "udscbt")
    return J*(t*t-c*c)*(t*t-u*u)*(c*c-u*u)*(b*b-s*s)*(b*b-d*d)*(s*s-d*d)/T**12
eta = 6.1e-10; nBs = eta/7.04
out=[]; P=lambda s="": (print(s), out.append(s))
P("Supply D(T) vs need:   need eta = 6.1e-10  (or n_B/s = 8.7e-11, a different bookkeeping)")
P(f"{'T (GeV)':>8s} {'D, PDG masses':>14s} {'D, running':>12s} {'gap (eta/D)':>12s}")
for T in (160, 131, 100, 50, 25, 18, 10):
    P(f"{T:8.0f} {D(T,pole):14.1e} {D(T,run):12.1e} {eta/D(T,pole):12.1e}")
P("\nHow many powers of ten each premise controls (headline T = 100 GeV, PDG masses: gap = %.1e = 10^%.1f)"
  % (eta/D(100,pole), np.log10(eta/D(100,pole))))
for lab,T in (("P1 creation temperature 131 GeV (sphaleron switch-off, lattice)",131),("P1 ... 160 GeV (crossover)",160)):
    P(f"  {lab:58s}: gap 10^{np.log10(eta/D(T,pole)):.1f}")
P(f"  {'P2 masses at creation energy instead of PDG':58s}: gap 10^{np.log10(eta/D(100,run)):.1f}  (larger)")
P(f"  {'need written as n_B/s instead of eta':58s}: gap 10^{np.log10(nBs/D(100,pole)):.1f}")
for tgt,lab in ((eta,"eta"),(nBs,"n_B/s")):
    Tc = brentq(lambda T: D(T,pole)-tgt, 1, 200)
    P(f"  temperature at which the measured bias would be enough ({lab}): T = {Tc:.0f} GeV")
P("  P3 the (m/T)^12 form: valid only when all masses << T; at T ~ 20 GeV top (173) and bottom (4) are not << T ->")
P("     the estimate itself breaks where it would 'succeed' (the real suppression is NOT this simple power law there)")
P("  P4 zero starting excess: if false, the needed creation is 0 -> the whole gap (10^9.7) disappears")
P("  P5 out-of-balance: smooth crossover (lattice, Higgs 125 GeV) multiplies supply by ~0 -> independent of P1-P3")
open("RESULT_127_numbers.txt","w").write("\n".join(out)+"\n")
