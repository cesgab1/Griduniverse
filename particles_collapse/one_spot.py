"""
Why each particle lands in one spot: every candidate grid rule tested against data, then the Born rule.
Part 1  Commit rules (the growing grid picks ONE branch when the branches differ by 'one unit'):
  A  spacetime volume   : dV >= 1 grid cell                     (Gm/c^2) d^2 c tau / l_d^4 >= 1
  B  clock mismatch     : proper-time difference >= 1 grid tick  E_G tau >= m c l_d
  C  action difference  : gravitational action >= hbar, mass resolved at the grid scale   (= Diosi-Penrose with R0 = l_d)
  D  action difference  : as C, but mass smeared over R0 >= 0.5 Angstrom (the allowed Diosi-Penrose)
Part 2  The odds: why |psi|^2?
  (a) Frequencies: weight of 'wrong-frequency' branches under |psi|^2 vs under plain branch counting
  (b) Triple slit: a path-sum with measure |psi|^n gives third-order interference kappa = 0 only for n = 2; the Sinha 2010 bound
      |kappa| < 1e-2 (relative to two-path interference) -> measured exponent n
"""
import numpy as np
from scipy.stats import binom
from scipy.optimize import brentq
G, hbar, c, amu = 6.6743e-11, 1.054571817e-34, 2.99792458e8, 1.66054e-27
ld = 0.603 * 1.616255e-35
def EG_sphere(m, R, d):
    lam = d / R; W = G*m*m/d if lam >= 2 else G*m*m/R*(6/5 - lam**2/2 + 3/16*lam**3 - lam**5/160)
    return 6/5*G*m*m/R - W

print("=== Part 1: can the grid pick the spot? ===")
# data: electron double slit (1e-6 m, 1e-8 s); Kovachy 2015 Rb-87 atoms 0.54 m apart for ~1 s; C60 (1e-7 m, 1e-2 s)
m_e, m_rb = 9.109e-31, 87*amu
tests = [("electron double slit", m_e, 1e-6, 1e-8, 2.8e-15), ("Rb atom, 54 cm for 1 s (Kovachy 2015)", m_rb, 0.54, 1.0, 5e-15),
         ("C60 molecule", 720*amu, 1e-7, 1e-2, 5e-15)]
print("rule A (one cell of spacetime volume): collapse after")
for n, m, d, t, R in tests:
    t_c = ld**4 / ((G*m/c**2) * d*d * c)       # time at which dV reaches one cell
    print(f"   {n:40s} {t_c:9.1e} s   (superposition lasted {t:.0e} s) -> {'RULED OUT' if t_c < t else 'ok'}")
print("rule B (proper-time mismatch of one grid tick): collapse after")
for n, m, d, t, R in tests:
    E = EG_sphere(m, R, d) if m > 1e-28 else G*m*m/R      # nucleon-scale mass concentration
    t_c = m*c*ld / E
    print(f"   {n:40s} {t_c:9.1e} s   (lasted {t:.0e} s) -> {'RULED OUT' if t_c < t else 'ok'}")
print("rule C (action difference = hbar, mass resolved to the grid spacing):")
R0_bound = 0.54e-10
print(f"   collapse-induced radiation scales as R0^-3; Gran Sasso (Donadi et al. 2021) requires R0 > ~{R0_bound:.1e} m;"
      f" grid-scale R0 = {ld:.1e} m overshoots the measured limit by ~10^{3*np.log10(R0_bound/ld):.0f} -> RULED OUT")
print("rule D (action difference = hbar, R0 >= 0.5 Angstrom): allowed today, but R0 is not a grid quantity.")
for n, m, d, t, R in tests:
    E = EG_sphere(m, max(R, R0_bound), d)
    print(f"   {n:40s} collapse after {hbar/E:9.1e} s -> {'RULED OUT' if hbar/E < t else 'ok'}")
m_bmv = 1e-14; R_bmv = (3*m_bmv/(4*np.pi*3500))**(1/3)
print(f"   decisive case: 1e-14 kg microdiamond, 250 um apart -> collapse after {hbar/EG_sphere(m_bmv, R_bmv, 250e-6):.3f} s (BMV holds it 2.5 s)")

print("\n=== Part 2(a): why |psi|^2 and not branch counting ===")
p, eps = 0.9, 0.05          # e.g. a detector pixel with 90% of the intensity, repeated N times
for N in (10, 100, 1000, 10000):
    k = np.arange(N + 1); near = np.abs(k/N - p) <= eps
    w_born = binom.pmf(k, N, p)[near].sum()
    logC = np.array([np.sum(np.log(np.arange(1, N+1))) - np.sum(np.log(np.arange(1, kk+1))) - np.sum(np.log(np.arange(1, N-kk+1))) for kk in k])
    frac_count = np.exp(logC[near] - N*np.log(2)).sum()
    print(f"   N = {N:6d}: share of the state (|psi|^2 weight) with observed frequency 0.90 ± 0.05 = {w_born:.6f};  share of branches by counting = {frac_count:.1e}")
print("   -> with |psi|^2 weights, almost all of the state records the 90% frequency that experiments see; counting branches equally")
print("      predicts 50% every time and is contradicted by any unequal-intensity experiment. The grid's path-sum must use |psi|^2.")

print("\n=== Part 2(b): the triple slit measures the exponent in the Born rule ===")
rng = np.random.default_rng(3)
def kappa_ratio(nexp, trials=4000):
    out = []
    for _ in range(trials):
        a = rng.uniform(0.6, 1.0, 3) * np.exp(1j*rng.uniform(0, 2*np.pi, 3))
        P = lambda *s: abs(sum(a[i] for i in s))**nexp
        k = P(0, 1, 2) - P(0, 1) - P(0, 2) - P(1, 2) + P(0) + P(1) + P(2)
        two = abs(P(0, 1) - P(0) - P(1))
        out.append(abs(k) / max(two, 1e-9))
    return np.median(out)
for nexp in (2.0, 2.02, 2.1, 3.0, 4.0):
    print(f"   measure |psi|^{nexp:<4}: third-order interference / two-path interference = {kappa_ratio(nexp):.2e}")
nmax = brentq(lambda x: kappa_ratio(2 + x) - 1e-2, 1e-4, 0.5)
print(f"   Sinha et al. 2010: ratio < 1e-2  ->  Born-rule exponent n = 2 within ± {nmax:.3f}")
