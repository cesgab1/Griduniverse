"""
ITERATION 2: where do the kicks come from, quantum mechanically?

Take the grid's link field phi as a free massless QUANTUM field (the long-wavelength limit of the quantised tension
links; cf. ../area_law.py). A domain's accumulated tension is the field averaged over the domain (Gaussian window, radius
R) and accumulated with memory rate gamma:
    tau(T) = ∫_{-inf}^{T} e^{-gamma (T - t)} phi_R(t) dt
Its variance (symmetrised quantum correlator), per mode k:
    Var = ∫ d^3k/(2pi)^3 · [coth(k / 2 Theta) / (2k)] · e^{-k^2 R^2} · 1/(gamma^2 + k^2)
          (vacuum: coth -> 1; relic thermal background of the field's own quanta at temperature Theta: full coth)
Question 1: does the vacuum alone give a random walk (Var ∝ 1/gamma, i.e. ∝ memory time)? Or only a log?
Question 2: does a thermal background give Var ∝ Theta/gamma, independent of domain size R?
Expanding universe: for a CONFORMAL field (photon-like; conformally coupled scalar xi = 1/6), phi = chi/a, dt = a d(eta); chi lives in flat space
at the fixed COMOVING temperature Theta_c, and memory gamma = kappa H in time is kappa * a * H = kappa * adot in eta.
=> Var_thermal = Theta_c / (4 pi kappa adot)  ∝ 1/adot (adot = da/dt): the counting of Claim 1 with NO comoving length put in.
Validity (independent check, Oct 2026): (a) the adiabatic formula needs kappa >> |q| (else use the moment ODE of iteration 1:
d<tau^2>/dt = -2 kappa H <tau^2> + Theta_c/(2 pi a)); (b) white noise needs kappa H << Theta (true for the CMB);
(c) NOT valid for gravitons or minimally coupled fields (they feel a''/a at horizon scale, where the relevant modes sit);
(d) a conformal field keeps its Planck spectrum with Theta ∝ 1/a (no particle creation). Photons are a vector: scalar = analogy.
Units: hbar = c = 1; numbers below in units of 1/R.
"""
import numpy as np
from scipy.integrate import quad
def var(gamma, Theta, R=1.0, part="total"):
    def f(k):
        w = np.exp(-k*k*R*R)/(gamma*gamma + k*k)*k/(4*np.pi**2)          # 4 pi k^2 /(2pi)^3 /(2k) = k/(4 pi^2)
        if part == "vac": return w
        x = k/(2*Theta); occ = 1/np.tanh(x) if x > 1e-8 else 1/x
        return w*(occ - 1) if part == "therm" else w*occ
    pts = sorted({gamma, 1/R, 2*Theta}); s = 0; lo = 0
    for p in pts + [40/R]:
        if p > lo: s += quad(f, lo, p, limit=400, epsabs=0, epsrel=1e-10)[0]; lo = p
    return s
out = ["ITERATION 2: quantum origin of the kicks (massless quantum field, domain radius R = 1, memory rate gamma)", ""]
out.append("Q1  VACUUM only: accumulated variance vs memory time 1/gamma")
gs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
v = np.array([var(g, 1.0, part="vac") for g in gs])
for g, x in zip(gs, v): out.append(f"   memory 1/gamma = {1/g:8.0e}:  Var = {x:.5f}    exact e^(g^2R^2)E1(g^2R^2)/(8pi^2) = {np.exp(g*g)*__import__('scipy.special',fromlist=['exp1']).exp1(g*g)/(8*np.pi**2):.5f}")
ex = np.diff(np.log(v))/np.diff(np.log(1/gs)); out.append(f"   growth exponent d lnVar / d ln(memory): {', '.join(f'{e:.3f}' for e in ex)}   (random walk would be 1)")
out.append("")
out.append("Q2  THERMAL background (thermal part only): Var vs memory, temperature and domain size")
rows = []
for Th in (1e-2, 1e-1, 1.0):
    for g in (1e-3, 1e-4, 1e-5):
        for R in (0.5, 1.0, 4.0):
            x = var(g, Th, R=R, part="therm"); rows.append((Th, g, R, x))
            out.append(f"   Theta {Th:5.0e}  memory {1/g:6.0e}  R {R:3.1f}:  Var = {x:11.4e}   Theta/(4 pi gamma) = {Th/(4*np.pi*g):11.4e}   ratio {x/(Th/(4*np.pi*g)):.4f}")
out.append("   -> Var = Theta/(4 pi gamma) to the accuracy shown once memory >> R and >> 1/Theta: independent of domain size R")
out.append("")
out.append("Q3  Expanding universe (conformal field): Var = Theta_c/(4 pi kappa adot); kick rate in physical time = Theta(t)/(2 pi), Theta ∝ 1/a")
out.append("   = iteration 1's kick rate ∝ 1/a, now DERIVED; the 'comoving length' of law_from_grid.md section 7 is not needed:")
out.append("   the 1/a comes from the redshifting temperature of real relic quanta of a CONFORMAL field (photons; not gravitons).")
c = 2.998e8; H0 = 67.5e3/3.0857e22; kB, hb = 1.380649e-23, 1.0546e-34
for lab, T in (("CMB photons 2.725 K (kappa = 1, window R = Planck length)", 2.725),):
    Th = kB*T/hb                                       # temperature as an angular frequency
    out.append(f"   {lab}: Theta/H0 = {Th/H0:.1e} -> thermal Var/vacuum Var ≈ {Th/(4*np.pi*H0)/(np.log(c/(H0*1.6e-35))/(4*np.pi**2)):.1e}")
out.append("Reading: the VACUUM alone gives only a logarithm (fluctuations of a massless field are anti-correlated in time: no")
out.append("random walk; would give a near-constant dark energy). REAL thermal quanta give an honest random walk whose size ∝ 1/adot.")
txt = "\n".join(out); print(txt); open("iter2_quantum_kicks.txt", "w").write(txt + "\n")
