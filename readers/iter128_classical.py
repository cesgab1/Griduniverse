"""PREREG 128: classical tools (Newton shell theorem + heat physics) applied to matter creation at 16-20 GeV."""
import numpy as np
G=6.674e-11; hbar=1.0546e-34; c=2.998e8; GeV=1.602e-10   # J
out=[]; P=lambda s="": (print(s), out.append(s))
def gstar(T):   # counted relativistic species (bosons + 7/8 fermions) from measured particle list
    b = 2+16; f = 3*4 + 3*2 + 5*12            # photon, gluons; e mu tau; 3 nu; u d s c b quarks
    if T > 80: b += 6+3+1; f += 12            # W, Z, Higgs, top present at ~100+ GeV (crude step)
    return b + 7/8*f
def H(T):       # Newton: H^2 = 8 pi G rho/3, rho = energy density / c^2
    u = np.pi**2/30*gstar(T)*(T*GeV)**4/(hbar*c)**3     # J/m^3
    return np.sqrt(8*np.pi*G*(u/c**2)/3)
P("C1 timeline (radiation era, t = 1/(2H)):")
for T in (160,100,20,16):
    P(f"   T = {T:4d} GeV ({T*GeV/1.380649e-23:.1e} K): g* = {gstar(T):6.2f}, time = {1/(2*H(T)):.1e} s")
v, g, B = 246.0, 0.65, 1.9
E = 4*np.pi*v*B/g
P(f"\nC2 hopping: barrier E = 4 pi v B / g = {E/1000:.1f} TeV  (B range 1.5-2.7 -> {4*np.pi*v*1.5/g/1000:.1f}-{4*np.pi*v*2.7/g/1000:.1f} TeV)")
for T in (20,16):
    lnr = np.log(T*GeV/(hbar*H(T)))
    P(f"   T = {T} GeV: barrier/T = {E/T:.0f}; hops keep up only if barrier/T < ln(T/hbar H) = {lnr:.1f} (+/-10)")
    P(f"      -> hop rate / expansion ~ 10^{(lnr - E/T)/np.log(10):.0f};  largest allowed barrier = {lnr*T/1000:.2f} TeV"
      f" ({(lnr-10)*T/1000:.2f}-{(lnr+10)*T/1000:.2f}); needs barrier {E/(lnr*T):.0f}x lower than measured inputs give")
aw = g**2/(4*np.pi)
P("\nC3 balance at 16-20 GeV (weak rate ~ alpha_w^2 T vs expansion):")
for T in (20,16):
    P(f"   T = {T} GeV: weak rate / expansion = {aw**2*T*GeV/(hbar*H(T)):.1e}")
P("   particles that stop being produced between 16 and 20 GeV: none (W, Z, top, Higgs already gone at >~80 GeV;")
P("   bottom 4.2 GeV, charm 1.3, tau 1.8 still fully in balance) -> nothing out of balance in the measured particle list")
open("RESULT_128_numbers.txt","w").write("\n".join(out)+"\n")
