"""PREREG 126: anti-readers. (a) capacity equality, (b) CP asymmetry vs mass, (c) matter-excess budget."""
import numpy as np
from scipy.stats import spearmanr
out=[]; P=lambda s="": (print(s), out.append(s))
P("(a) particle vs antiparticle capacity (mass) equality:")
for lab,v in (("neutral kaon mass |dm|/m (90%)",6e-19),("proton/antiproton q/m (BASE 2022)",1.6e-11),
              ("electron/positron g",2.1e-12),("electron/positron mass (90%)",8e-9)):
    P(f"    {lab:36s} equal to {v:.0e}")
sysm = {"K0":0.497611,"D0":1.86484,"B0":5.27972,"Bs":5.36691,"Lambda_b":5.61957}
head = {"K0":2.228e-3,"D0":15.4e-4,"B0":0.691,"Bs":0.236,"Lambda_b":0.0245}     # |eps_K|, |dA_CP|, sin2beta, A_CP(Kpi), A_CP
alt  = dict(head, K0=1.67e-3*2.228e-3, B0=0.0824)                                  # eps'/eps*|eps| ~ direct; B0 direct Kpi
P("\n(b) CP asymmetry vs mass (capacity):")
for lab,d in (("headline (largest well-measured per system)",head),("alternative (direct asymmetries)",alt)):
    keys=list(sysm); rho,p = spearmanr([sysm[k] for k in keys],[abs(d[k]) for k in keys])
    P(f"    {lab:44s}: Spearman rho = {rho:+.2f}, p = {p:.2f} -> {'TREND' if p<0.05 else 'no significant trend'}")
P("    power note: with 5 systems only rho = +/-1.0 reaches p < 0.05; three systems sit at 5.3-5.6 GeV (all b-quark)")
mu,md,ms,mc,mb,mt = 0.00216,0.00470,0.0935,1.2730,4.183,172.57
J=3.18e-5; T=100.0
D = J*(mt**2-mc**2)*(mt**2-mu**2)*(mc**2-mu**2)*(mb**2-ms**2)*(mb**2-md**2)*(ms**2-md**2)/T**12
eta=6.1e-10
P(f"\n(c) matter excess: standard CP measure J*prod(dm^2)/T^12 at T=100 GeV = {D:.1e}; needed eta = {eta:.1e};"
  f" short by {eta/D:.0e}")
open("RESULT_126_numbers.txt","w").write("\n".join(out)+"\n")
