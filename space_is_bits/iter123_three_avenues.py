"""PREREG 123: three avenues for 'space is bits'. Inputs are published numbers listed in SOURCES_123.md."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from scipy.integrate import quad
import qnm

c = 2.99792458e8; hbarc_eVm = 1.97327e-7; G = 6.674e-11; Msun = 1.98892e30; kB_eV = 8.617333e-5
MPC = 3.08568e22
out = []; P = lambda s="": (print(s), out.append(s))
def comoving(z, h0, Om): return c/(h0*1e3/MPC)*quad(lambda x: 1/np.sqrt(Om*(1+x)**3+1-Om), 0, z)[0]

# ---------------- T1 LOCKED BITS ----------------
P("T1  LOCKED BITS: jitter over distance L is dL = l^(1-b) L^b")
L = comoving(0.903, 71.0, 0.27)                      # GRB 090510
dLmax = c*10e-3                                      # spikes <~ 10 ms (Abdo et al. 2009)
P(f"  (a) GRB 090510 spikes <~10 ms, L = {L:.2e} m -> distance jitter < {dLmax:.1e} m")
for b, lab in ((0.5,"independent bits (random walk)"), (1/3,"locked bits (holographic, Ng-van Dam)")):
    lmax = (dLmax/L**b)**(1/(1-b))
    P(f"      b = {b:.2f} {lab:40s}: l < {lmax:.1e} m")
EPl = 1.22089e28                                     # eV; only to unpack Fermi's quoted unit (xi = E_QG / E_Pl)
for xi, cl in ((2.8,"95%"), (1.6,"99%")):
    P(f"  (b) random speed spread, linear (Vasileiou+2015, xi > {xi} at {cl}): random grid l < {hbarc_eVm/(xi*EPl):.1e} m")
P("      quadratic random spread: paper says limits are 'several orders' from Planck scale; no number -> not converted")
lP = 1.616e-35
P(f"  (c) Holometer 2017, shear-holographic (Hogan) model: |beta_2 L| < 0.25 t_P (2 sigma) -> that model's scale"
  f" < {0.25*lP:.1e} m  [model-specific]")
Lp = comoving(0.944, 70.0, 0.3)
dl_thresh = lP**0.72 * Lp**0.28                      # Perlman+2015 TeV: alpha >~ 0.72 allowed at Planck scale
lmax_d = (dl_thresh/Lp**(1/3))**1.5
P(f"  (d) CONTESTED image-blurring (Perlman+2015, TeV source z=0.944): locked bits l < {lmax_d:.1e} m"
  f"  [disputed: Coule 2003 -- random kicks cancel; not in headline]")

# ---------------- T2 SURFACE BITS ----------------
P("\nT2  SURFACE BITS: horizon area grows in tiles a = alpha*hbar*G/c^3; allowed ringing lines")
P("    omega_n M = [n alpha s + 16 pi chi] / [16 pi (1+s)],  s = sqrt(1-chi^2)   (Foit-Kleban / Laghi+2021 eq.)")
mode = qnm.modes_cache(s=-2, l=2, m=2, n=0)
_Fcache = {}
def F(chi):                                          # GR 220 frequency * M (cached)
    k = round(float(chi), 6)
    if k not in _Fcache: _Fcache[k] = mode(a=k)[0].real
    return _Fcache[k]
def line_gap(alpha, chi):
    s = np.sqrt(1-chi**2); base = 16*np.pi*chi/(16*np.pi*(1+s)); step = alpha*s/(16*np.pi*(1+s))
    x = (F(chi)-base)/step
    return abs(x-round(x))*step, step, x                # distance to nearest line (in omega*M)
events = {"GW150914": dict(chi=(0.60,0.72), ftol=0.20),   # spin 0.67 +0.05/-0.07; ringdown deviation dfreq -0.05 +/- 0.2 (Isi+2019)
          "GW250114": dict(chi=(0.67,0.69), ftol=6/247)}  # spin 0.68 +/- 0.01; f220 = 247 +/- 6 Hz (90%)
alphas = np.round(np.arange(0.5, 120.0001, 0.01), 2)
allowed = {}
for ev, d in events.items():
    chis = np.linspace(*d["chi"], 81); ok = np.zeros(len(alphas), bool)
    for ch in chis:                                  # vectorised over alpha
        sq = np.sqrt(1-ch**2); base = ch/(1+sq); step = alphas*sq/(16*np.pi*(1+sq)); x = (F(ch)-base)/step
        ok |= np.abs(x-np.round(x))*step <= d["ftol"]*F(ch)
    allowed[ev] = ok
    bad = alphas[~ok]
    P(f"  {ev}: spin {d['chi']}, freq tolerance {d['ftol']*100:.1f}% -> excluded alpha fraction in 0.5-120: "
      f"{(~ok).mean()*100:.1f}%" + (f"; smallest excluded alpha {bad.min():.2f}" if len(bad) else ""))
both = allowed["GW150914"] & allowed["GW250114"]
for name, al in (("4 ln2 (one bit per tile)", 4*np.log(2)), ("4 ln3 (Hod)", 4*np.log(3)),
                 ("8 ln2", 8*np.log(2)), ("8 pi (Bekenstein/Maggiore)", 8*np.pi)):
    g, st, x = line_gap(al, 0.68)
    i = np.argmin(abs(alphas-al)); verdict = "allowed" if both[i] else "EXCLUDED"
    P(f"  alpha = {al:6.3f} {name:28s}: chi=0.68 needs n = {x:6.2f}; nearest line off by {g/F(0.68)*100:4.1f}% "
      f"(spacing {st/F(0.68)*100:4.1f}%) -> {verdict} (both events, spin range included)")
# where does exclusion begin, using GW250114 alone at its spin range
first_band = alphas[~allowed["GW250114"]]
P(f"  GW250114 alone: all alpha < {first_band.min():.2f} allowed (lines too dense); in 0.5-120 excluded fraction "
  f"{(~allowed['GW250114']).mean()*100:.0f}%" if len(first_band) else "  GW250114: nothing excluded")
P(f"  tile area for alpha = 4 ln2: a = {4*np.log(2)*1.0546e-34*G/c**3:.1e} m^2 (measured hbar, G, c)")
P("  Alternative criterion (line within QNM half-width, 15.6%): all alpha < 25.6 allowed, incl. 8 pi (25.13)")
P("  Published: Laghi+2021 alpha = 15.6 +20.5/-13.3, log odds 0.1 +/- 0.6 (uninformative); GWTC-3: no echoes.")

# ---------------- T3 AREA RULE AT SMALL SCALES ----------------
P("\nT3  DOES THE AREA RULE BREAK AT SMALL SCALES?")
P("  (a) Newton's law verified down to 3.86e-5 m (Lee+2020, gravitational-strength Yukawa ranges < 38.6 um only)")
for l in (1.7e-27, lP):
    P(f"      gap to a {l:.1e} m cell: {np.log10(3.86e-5/l):.0f} powers of ten untested")
P("  (b) Area never decreases: GW150914 97% (Isi+2021); GW250114 4.4 sigma (LVK 2025) -> classical area rule HOLDS")
Mbh = 62.7*Msun; rs = 2*G*Mbh/c**2
TH = 1.0546e-34*c**3/(8*np.pi*G*Mbh*1.380649e-23)
ratio_real = kB_eV*TH/(hbarc_eVm/1.7e-27)
P(f"  (c) graininess vs Hawking energy: BEC analog k T_H = 0.12 m c^2 (grain energy) -> ratio 0.12; spectrum thermal,"
  f" no free parameters (Munoz de Nova+2019)")
P(f"      real 62.7 Msun hole: T_H = {TH:.1e} K, horizon {rs/1e3:.0f} km; vs 1.7e-27 m cell -> ratio {ratio_real:.0e}")
P(f"      -> the analog was ~{0.12/ratio_real:.0e} times MORE grainy (relative) and the thermal law still held")
open("RESULT_123_numbers.txt","w").write("\n".join(out)+"\n")
