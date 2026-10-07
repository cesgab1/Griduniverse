"""PREREG 125: does capacity enter quantum reading only through h/m?  eps = 2 * d(alpha)/alpha."""
import numpy as np
ref = {"Fan 2023 g-2 (headline)": (137.035999166, 15e-9), "Hanneke 2008 g-2 (alt)": (137.035999084, 51e-9)}
atoms = {"Cs-133 (Parker 2018)": (137.035999046, 27e-9, 132.905), "Rb-87 (Morel 2020)": (137.035999206, 11e-9, 86.909)}
out=[]; P=lambda s="": (print(s), out.append(s))
for rn,(r,sr) in ref.items():
    P(f"reference: {rn}  alpha^-1 = {r:.9f}({sr*1e9:.0f})")
    signs=[]; sig=[]
    for an,(a,sa,m) in atoms.items():
        dalpha = -(a-r)/a                      # relative change in alpha
        s = abs(a-r)/np.hypot(sa,sr); eps = 2*dalpha
        signs.append(np.sign(eps)); sig.append(s)
        P(f"   {an:22s} mass {m:7.3f} u: capacity correction eps(h/m) = {eps:+.2e}  ({s:.1f} sigma)")
    opp = signs[0]!=signs[1]; both2 = min(sig)>2
    verdict = ("single monotonic capacity law EXCLUDED (opposite signs, each > 2 sigma)" if opp and both2 else
               "inconclusive (opposite signs but not both > 2 sigma)" if opp else "same sign -> capacity law allowed")
    P(f"   -> {verdict}")
d = (137.035999206-137.035999046)/np.hypot(11e-9,27e-9)
P(f"\nCs vs Rb directly: {d:.1f} sigma (papers: '> 5 sigma')")
P("Neutron h/m (Kruger 1999) relative uncertainty 7.3e-8: ~100x too coarse to see these 1e-9 effects.")
open("RESULT_125_numbers.txt","w").write("\n".join(out)+"\n")
