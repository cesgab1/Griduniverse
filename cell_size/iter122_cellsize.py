"""PREREG 122: grid cell size from GRB photon timing. No Planck length used in any calculation."""
import numpy as np
from scipy.integrate import quad

hbarc = 1.97327e-16            # GeV*m (measured)
MPC = 3.08568e22               # m
def H0s(h0): return h0*1e3/MPC # 1/s
def K(n, z, Om, OL): return quad(lambda x: (1+x)**n/np.sqrt(Om*(1+x)**3+OL), 0, z)[0]
def EQG_from_dt(n, Eh, El, dt, z, h0, Om, OL):
    return ((1+n)/2*(Eh**n-El**n)*K(n,z,Om,OL)/H0s(h0)/dt)**(1/n)
def dt_from_EQG(n, Eh, El, EQG, z, h0, Om, OL):
    return (1+n)/2*(Eh**n-El**n)/EQG**n*K(n,z,Om,OL)/H0s(h0)
l_lin  = lambda E1: hbarc/E1
l_quad = lambda E2, f: np.sqrt(12/f)*hbarc/E2
l_to_E2 = lambda l, f: np.sqrt(12/f)*hbarc/l

out = []
P = lambda s="": (print(s), out.append(s))

# ---- A. GRB 090510 (Fermi) ----
zA, cosA = 0.900, (71.0, 0.27, 0.73)
EhA = 28.0
rows = {"(a) any <1 MeV":(859e-3,1e-4), "(b) main <1 MeV":(299e-3,1e-4),
        "(c) main >0.1 GeV":(181e-3,0.1), "(d) >1 GeV":(99e-3,1.0)}
E1a = EQG_from_dt(1, EhA, 1e-4, 0.859, zA, *cosA)
EPl_ref = 1.22089e19   # used ONLY to compare with the paper's printed ratio (C1)
ratio = E1a/EPl_ref
P(f"C1 reproduction: row (a) -> E_QG,1 = {E1a:.3e} GeV = {ratio:.3f} x paper unit (paper 1.19); "
  f"{'PASS' if abs(ratio/1.19-1)<0.05 else 'FAIL'}")
assert abs(ratio/1.19-1) < 0.05

P("\nA. GRB 090510 (Fermi, single 28 GeV photon, z = 0.900)")
P(f"{'start assumption':22s} {'dt limit':>9s} {'linear l <':>11s} {'grid(diag) l <':>15s} {'grid(axis) l <':>15s}")
for k,(dt,El) in rows.items():
    E1 = EQG_from_dt(1, EhA, El, dt, zA, *cosA); E2 = EQG_from_dt(2, EhA, El, dt, zA, *cosA)
    P(f"{k:22s} {dt*1e3:7.0f}ms {l_lin(E1):11.2e} {l_quad(E2,1/3):15.2e} {l_quad(E2,1):15.2e}")
E1g = K(1,zA,*cosA[1:])/H0s(cosA[0])/0.030
P(f"{'(g) spike lags':22s} {'30ms/GeV':>9s} {l_lin(E1g):11.2e} {'--':>15s} {'--':>15s}")

# ---- B. GRB 221009A (LHAASO) ----
zB, cosB = 0.151, (67.36, 0.315, 0.685)
E1B, E2B = 1.0e20, 6.9e11
P("\nB. GRB 221009A (LHAASO, 0.2-7 TeV, z = 0.151, 95% ML subluminal)")
P(f"   linear l < {l_lin(E1B):.2e} m;  regular grid l < {l_quad(E2B,1/3):.2e} m (diagonal), {l_quad(E2B,1):.2e} m (axis)")
dtmaxB2 = dt_from_EQG(2, 7000, 200, E2B, zB, *cosB); dtmaxB1 = dt_from_EQG(1, 7000, 200, E1B, zB, *cosB)
P(f"   i.e. data allow at most {dtmaxB2:.1f} s (grid) / {dtmaxB1:.1f} s (linear) delay of 7 TeV behind 0.2 TeV photons")

# ---- E3: seconds for trial cell sizes (regular grid, diagonal, LHAASO burst) ----
P("\nE3. Delay of 7 TeV vs 0.2 TeV photons from GRB 221009A for a regular grid (diagonal):")
for l in (1e-24, 1e-26, 1e-27, 1e-28, 1e-30, 1e-35):
    P(f"   cell {l:.0e} m -> {dt_from_EQG(2,7000,200,l_to_E2(l,1/3),zB,*cosB):.3g} s")

# ---- free choice 3: cosmology spread ----
P("\nCosmology spread (H0 67-73, Om 0.27-0.32) on headline grid bound:")
vals = [l_quad(E2B,1/3)*np.sqrt(K(2,zB,Om,1-Om)/H0s(h)/(K(2,zB,0.315,0.685)/H0s(67.36))) for h in (67,73) for Om in (0.27,0.32)]
P(f"   l < {min(vals):.2e} to {max(vals):.2e} m  (bound scales as sqrt of the distance integral)")

lP = 1.616e-35
P(f"\nComparison only: Planck length {lP:.2e} m. Regular grid allowed down to {l_quad(E2B,1/3)/lP:.1e} x it; "
  f"linear graininess must be < {l_lin(E1B)/lP:.2f} x it.")
open("RESULT_122_numbers.txt","w").write("\n".join(out)+"\n")
