"""PREREG 130: classical opacity-limited fragmentation floor vs lowest observed masses."""
import numpy as np
G=6.674e-11; k=1.381e-23; mH=1.6735e-27; sigma=5.670e-8; Msun=1.989e30; MJ=1.898e27
out=[]; P=lambda s="": (print(s), out.append(s))
def mmin(T, mu=2.33, f=1.0):
    """heating of a collapsing Jeans lump (G M^2/R per free-fall time) = blackbody cooling f*4piR^2 sigma T^4,
    with R = G M mu mH/(kT). Gives M ~ T^(1/4) mu^(-9/4) G^(-3/2) f^(-1/2)  (Rees 1976 scaling)."""
    a = G*mu*mH/(k*T)
    return np.sqrt(G**1.5*np.sqrt(3/(4*np.pi))/(4*np.pi*sigma*f)) * a**(-9/4) / T**2
P("Classical floor (our derivation, order-unity factors):")
for T in (7,10,15,20):
    P(f"   T = {T:2d} K: {mmin(T)/MJ:.2f} Jupiter masses (f=1); {mmin(T,f=0.3)/MJ:.2f} (f=0.3, inefficient radiator)")
pub = 0.003*Msun/MJ
P(f"Check vs published: Whitworth 2018 0.003 +/- 0.001 Msun = {pub:.1f} +/- {0.001*Msun/MJ:.1f} MJup; "
  f"ours at 10 K = {mmin(10)/MJ:.1f} -> ratio {pub/(mmin(10)/MJ):.1f} (order-unity prefactors; PASS if within x3)")
floor_lo, floor_hi = 1.0, 4.0   # published range 0.001-0.004 Msun (Whitworth & Stamatellos 2006) -> ~1-4 MJup
obs = [("IC 348 (JWST NIRSpec, Luhman & Alves de Oliveira 2025)", 2.0, "spectroscopic members"),
       ("NGC 1333 (JWST NIRISS, Langeveld 2024)", 4.5, "none below ~4-5 MJup despite sensitivity"),
       ("NGC 2024 (JWST NIRCam, 2024)", 3.0, "photometric; turnover below ~12 MJup"),
       ("Upper Sco (Miret-Roig 2022)", 4.0, "photometric/kinematic"),
       ("Orion Trapezium JuMBOs (Pearson & McCaughrean 2023)", 0.6, "DISPUTED: NIRSpec follow-up -> background")]
P(f"\nPublished floor range ~{floor_lo:.0f}-{floor_hi:.0f} MJup; violation needs confirmed objects < floor/2 in >= 2 regions")
viol = 0
for lab, m, note in obs:
    below = m < floor_lo/2 if "DISPUTED" not in note else False
    viol += below
    P(f"   {lab:55s} lowest ~{m:.1f} MJup ({note}) -> {'below floor/2' if below else 'at/above floor'}")
P(f"Regions with confirmed objects below floor/2: {viol} -> {'VIOLATION' if viol>=2 else 'no violation'}")
open("RESULT_130_numbers.txt","w").write("\n".join(out)+"\n")

# ---- POST-HOC diagnosis (added after the pre-registered check FAILED; not part of the prereg) ----
P("\nPOST-HOC: the floor scales as (lump radius factor)^(-9/4), so the radius convention dominates.")
def mmin_beta(T, beta, mu=2.33, f=1.0):
    a = beta*G*mu*mH/(k*T)
    return np.sqrt(G**1.5*np.sqrt(3/(4*np.pi))/(4*np.pi*sigma*f)) * a**(-9/4) / T**2
for beta, lab in ((1.0,"crude R = GM mu mH/kT (pre-registered run)"),(0.4,"R = 2/5 GM mu mH/kT"),
                  (0.2,"virial Jeans lump R = GM mu mH/(5kT)")):
    P(f"   {lab:45s}: floor at 10 K = {mmin_beta(10,beta)/MJ:.1f} MJup")
