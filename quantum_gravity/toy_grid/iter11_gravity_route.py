"""
ITERATION 11 (Coalesce): the jostling reaches the tension THROUGH GRAVITY - like a rubber sheet stretched from both ends
that is pulled down, up and sideways, making waves via the pull.
Big advantage: gravity's strength is KNOWN (G), so there is no free coupling: the height becomes a prediction.
Order-of-magnitude estimate (units eV, hbar = c = 1):
  glow energy density rho_g <= Delta N_eff-limited value; its thermal fluctuations in a volume V: delta rho = sqrt(4 rho_g T / V)
  potential they make on scale L (Poisson): Phi_L ~ (3/2) delta rho / rho_crit x (L H)^2 ; largest scales dominate (L = 1/H)
  the expansion-rate route (delta K ~ dPhi/dt) gives no random walk (a time derivative, iteration 3); the time-stamp route
  (tick rate 1 + Phi) couples to Phi itself and is used here.
  most generous response of a Planck-strength tension to metric jitter: rho_DE ~ rho_crit x <Phi^2> x N_layers
"""
import numpy as np
H0 = 1.44e-33; rho_crit = 2.54e-11/0.69                     # eV, eV^4
Tnu = 1.945*8.617e-5; rho_nu1 = 7/8*np.pi**2/30*2*Tnu**4      # one neutrino-like species (2 dof x 7/8)
out = ["ITERATION 11: jostling through gravity (no free coupling)", ""]
for dN in (0.107, 0.01):
    rho_g = dN*rho_nu1/1.0; T = Tnu*(dN)**0.25                # glow energy and a representative temperature
    for N in (1e7, 1.2e15):
        rg, Tl = rho_g/N, T/N**0.25                             # per layer (each colder, Delta N_eff shared)
        dr = np.sqrt(4*rg*Tl*H0**3); Phi = 1.5*dr/rho_crit
        ratio = Phi**2*N/0.69
        out.append(f"Delta N_eff {dN:5.3f}, N = {N:7.1e} layers: Phi_H per layer ~ {Phi:.1e};  predicted rho_DE / observed ~ {ratio:.0e}")
out += ["",
        "-> Gravity is so weak, and the glow so faint, that gravity-mediated jostling makes a dark energy ~1e99-1e102 times too",
        "   small even with the most generous response. Smaller scales are worse (Phi_L ∝ sqrt(L)).",
        "   Gravitational WAVES instead: the relic-wave background is capped (Omega_GW <~ 1e-6, nucleosynthesis), and jostling",
        "   by gravity's own strength still brings a factor (Phi)^2 of order 1e-100.",
        "Lesson: the gravity route is PREDICTIVE (no free knob) and its prediction is far too small. The size of dark energy",
        "cannot come from jostling by the grid's own faint vibrations through any channel tried (direct, electrons, gravity).",
        "The jostling picture explains SHAPE, memory and smoothness; the SIZE must come from elsewhere - e.g. a grid tension that",
        "is already large and only MODULATED by the jostling (the jostling sets how it changes, not how big it is)."]
txt = "\n".join(out); print(txt); open("iter11_gravity_route.txt", "w").write(txt + "\n")
