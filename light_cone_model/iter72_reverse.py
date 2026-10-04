"""Iteration 72: run the Light-Cone Model backward from today's measured dark energy. See PREREG_72.md."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid as ctz
ODE, Ok, Or = 0.685, 0.0023, 9.1e-5
Om = 1 - ODE - Ok - Or
H0 = 67.4e3/3.0857e22            # 1/s
c, G, hbar = 2.998e8, 6.674e-11, 1.0546e-34
lP2 = hbar*G/c**3
Gyr = 3.156e16
rc0 = 3*H0**2/(8*np.pi*G)        # kg/m^3

def E_of(x, fading=True):        # x = ln a/a0
    z1 = np.exp(-x)
    base = Om*z1**3 + Or*z1**4 + Ok*z1**2
    if not fading: return np.sqrt(base + ODE)
    f = lambda lnE: np.exp(2*lnE) - base - ODE*(np.exp(lnE)*np.exp(x))**-0.5
    return np.exp(brentq(f, -80, 80))

x = np.linspace(-14, 40, 30001)
res = {}
for fading in (True, False):
    E = np.array([E_of(xi, fading) for xi in x])
    t = ctz(1/E, x, initial=0) + (1/E[0])/2*np.exp(-0) * 0   # starts ~0 at a=1e-6 (radiation era: t ~ 1/(2H))
    t += 1/(2*E[0])                                            # radiation-era offset
    t_s = t/H0
    rDE = ODE*(E*np.exp(x))**-0.5 if fading else ODE*np.ones_like(x)
    rM = Om*np.exp(-3*x)
    # event horizon (physical): a * int_t^inf c dt/a = (c/H0) e^x int_x^inf e^-x'/E dx'
    integ = np.exp(-x)/E
    tail = ctz(integ[::-1], -x[::-1], initial=0)[::-1]
    Reh = c/H0*np.exp(x)*tail
    RH = c/(H0*E)
    lnE = np.log(E); q = -1 - np.gradient(lnE, x)
    i0 = np.argmin(abs(x))
    ic = np.where(np.diff(np.sign(rDE - rM)))[0][0]
    ia = np.where(np.diff(np.sign(q)))[0][-1]
    res[fading] = dict(E=E, t=t_s, rDE=rDE, rM=rM, Reh=Reh, RH=RH, q=q, i0=i0, ic=ic, ia=ia)

out = ["Iteration 72: reverse calculation from today's measured dark energy (D0 = old capital C). PREREG_72.md", ""]
F = res[True]; L = res[False]
k = lambda rho_frac, R: rho_frac*rc0*R**2*G/c**2   # rho (as mass density) x R^2 in units c^4/G  -> rho c^2 R^2 G / c^4
for name, r in (("Fading (Claim 1)", F), ("Constant (Lambda)", L)):
    i0, ic, ia = r["i0"], r["ic"], r["ia"]
    out.append(f"{name}: age today {r['t'][i0]/Gyr:.2f} Gyr | crossing z = {np.exp(-x[ic])-1:.3f} at t = {r['t'][ic]/Gyr:.2f} Gyr "
               f"| acceleration begins z = {np.exp(-x[ia])-1:.3f} at t = {r['t'][ia]/Gyr:.2f} Gyr")
out.append("")
# R1
att = 3/(8*np.pi)*25/16
out.append(f"R1 matter x horizon^2 at the crossing (far-future quantum value {att:.4f} c^4/G):")
for name, r in (("Fading", F), ("Lambda", L)):
    ic = r["ic"]
    ke = k(r["rM"][ic], r["Reh"][ic]); kh = k(r["rM"][ic], r["RH"][ic])
    out.append(f"   {name}: event horizon {ke:.4f} ({ke/att:.2f} x) | Hubble radius {kh:.4f} ({kh/att:.2f} x)")
ke = k(F["rM"][F["ic"]], F["Reh"][F["ic"]]); kl = k(L["rM"][L["ic"]], L["Reh"][L["ic"]])
v1 = "PASS" if abs(ke/att-1) < 0.3 else "FAIL"
out.append(f"   Verdict (event horizon, fading): {v1}." + (" Lambda gives nearly the same -> no information." if abs(kl/ke-1) < 0.2 else ""))
# R2
zs = 1.86; xs = -np.log(1+zs); i_s = np.argmin(abs(x-xs))
out.append(f"R2 star-formation peak z = 1.86 at t = {F['t'][i_s]/Gyr:.2f} Gyr; crossing - peak = {(F['t'][F['ic']]-F['t'][i_s])/Gyr:.2f} Gyr; "
           f"acceleration start - peak = {(F['t'][F['ia']]-F['t'][i_s])/Gyr:.2f} Gyr (clue only)")
# R3
aHc0 = 1/np.sqrt(Ok)
D0 = ODE*np.sqrt(aHc0)                 # in units of critical density
rhoP = c**5/(hbar*G**2)                # Planck mass density kg/m^3 -> energy density rhoP c^2
D0_SI = D0*rc0*c**2                    # J/m^3
eV = 1.602e-19
scale = (D0_SI*(hbar*c)**3)**0.25/eV   # energy scale
rDE0_scale = (ODE*rc0*c**2*(hbar*c)**3)**0.25/eV
out.append(f"R3 D0 = {D0:.2f} x critical density = {D0_SI:.2e} J/m^3 = {D0*rc0/rhoP:.2e} of the Planck density; "
           f"energy scale D0^(1/4) = {scale*1e3:.2f} meV (today's dark energy: {rDE0_scale*1e3:.2f} meV).")
out.append("   Look-elsewhere flag: compare to our derived lightest neutrino 2.24 meV only as a curiosity, not evidence.")
# R4
aHc = F["E"]*np.exp(x)/F["E"][F["i0"]]*aHc0
out.append(f"R4 aH/c over history: min {aHc.min():.1f} (at z = {np.exp(-x[np.argmin(aHc)])-1:.2f}), today {aHc[F['i0']]:.1f}, "
           f"far future {aHc[-1]:.2e}. Empty light-cone state (aH/c = 1) reached: {'YES' if aHc.min() <= 1 else 'NEVER'}.")
# R5
ic = F["ic"]; R = F["Reh"][ic]
S_h = np.pi*R**2/lP2
s0 = 2.891e9                          # photon+neutrino entropy density today, k_B per m^3
Rcom = R*np.exp(-x[ic])
S_cmb = s0*4/3*np.pi*Rcom**3
S_bh = 3e104
out.append(f"R5 horizon entropy at crossing {S_h:.2e}; CMB+neutrino entropy inside it {S_cmb:.2e} (ratio {S_h/S_cmb:.1e}); "
           f"black holes today {S_bh:.0e} (ratio {S_h/S_bh:.1e}) -> " + ("PASS" if min(abs(np.log10(S_h/S_cmb)), abs(np.log10(S_h/S_bh))) < np.log10(3) else "FAIL"))
txt = "\n".join(out); print(txt); open("iter72_reverse.txt", "w").write(txt + "\n")

# consistency: today's rho_DE x R_eh^2 computed with the FADING history (what_sets_C.py used the Lambda history for this)
for name, r in (("Fading", F), ("Lambda", L)):
    i0 = r["i0"]; kt = k(r["rDE"][i0], r["Reh"][i0]); kc = k(r["rDE"][r["ic"]], r["Reh"][r["ic"]])
    out2 = f"Check {name}: dark energy x event horizon^2 today {kt:.4f} ({kt/att:.2f} of far-future), at crossing {kc:.4f}; " \
           f"event horizon today {r['Reh'][i0]/3.0857e25:.2f} Gpc"
    print(out2); open("iter72_reverse.txt", "a").write(out2 + "\n")
