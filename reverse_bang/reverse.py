"""Run the expansion backwards: size of today's visible universe, temperature, time -- classical vs grid density cap (bounce).
Calculation, not a data test. Radiation era with standard-model particle count g* (simplified steps); dark energy negligible early
(iteration 74: share ~3e-140 at Planck time)."""
import numpy as np
G, c, hb, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.3807e-23
yr, ly = 3.156e7, 9.461e15
R_obs = 46.5e9*ly                      # comoving radius of today's visible universe
T0 = 2.7255
rho_P = c**5/(hb*G**2); rho_cap = 0.41*rho_P   # grid / loop-quantum-cosmology style cap
lP = np.sqrt(hb*G/c**3)
def gstar(T):                           # effective particle count (rough steps)
    E = kB*T/1.602e-19                  # eV
    return 3.36 if E < 5e5 else 10.75 if E < 1e8 else 61.75 if E < 1e11 else 106.75
def gs(T): return 3.91 if kB*T/1.602e-19 < 5e5 else gstar(T)
def a_of_T(T): return T0/T*(3.91/gs(T))**(1/3)
def rho_of_T(T): return np.pi**2/30*gstar(T)*(kB*T)**4/(hb*c)**3/c**2
def t_of_T(T):                          # radiation era, classical: t = sqrt(3/(32 pi G rho))
    return np.sqrt(3/(32*np.pi*G*rho_of_T(T)))
rows = [("today", None, 13.8e9*yr, 1.0), ("fog clears (380,000 yr)", 3000, 3.8e5*yr, 1/1090), ("matter = radiation (~50,000 yr)", 9000, 5e4*yr, 1/3400)]
for lab, T in (("nuclei form (~3 min)", 1e9), ("neutrinos decouple (~1 s)", 1e10), ("collider limit (~1e-12 s)", 1.2e15),
               ("10^-30 s", None), ("density reaches the cap", "cap")):
    if T is None:
        # find T with t = 1e-30 s
        Ts = np.logspace(15, 33, 4000); tt = np.array([t_of_T(x) for x in Ts]); T = Ts[np.argmin(abs(np.log(tt/1e-30)))]
    if T == "cap":
        Ts = np.logspace(25, 33, 8000); r = np.array([rho_of_T(x) for x in Ts]); T = Ts[np.argmin(abs(np.log(r/rho_cap)))]
    rows.append((lab, T, t_of_T(T), a_of_T(T)))
print(f"{'stage':34s} {'temperature (K)':>16s} {'time':>12s} {'scale':>10s} {'visible-universe diameter today-patch':>40s}")
for lab, T, t, a in rows:
    d = 2*R_obs*a
    ds = f"{d/ly:.2e} light-years" if d > 1e13 else (f"{d/1.496e11:.2e} AU" if d > 1e9 else f"{d*1e3:.3f} mm" if d > 1e-6 else f"{d:.2e} m")
    ts = f"{t/yr:.2e} yr" if t > yr else f"{t:.1e} s"
    print(f"{lab:34s} {('-' if T is None else f'{T:.2e}'):>16s} {ts:>12s} {a:10.2e} {ds:>40s}")
Tcap = rows[-1][1]; a_cap = a_of_T(Tcap); D = 2*R_obs*a_cap
t_cl = t_of_T(Tcap)
# horizon at the cap (classical): ~ 2 c t ; number of disconnected regions in the patch
hor = 2*c*t_cl
print(f"\nAt the cap: density {rho_cap:.2e} kg/m^3 (0.41 Planck), temperature {Tcap:.1e} K, our visible universe was {D*1e3:.3f} mm across.")
print(f"Classical equation: density -> infinity a further {t_cl:.1e} s earlier (singularity). Grid cap: expansion rate falls to zero")
print(f"  there instead -- a BOUNCE (H^2 = 8piG/3 rho (1 - rho/rho_cap)): before it, a contracting phase, mirror image.")
print(f"Light-travel region at that moment ~{hor:.1e} m = {hor/lP:.1f} Planck lengths; regions in the patch that never touched: {(D/hor)**3:.1e}")
