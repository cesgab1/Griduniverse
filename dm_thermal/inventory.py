"""Infer what dark matter could be made of from what existed in the first seconds (known physics only).
Requirements: stable 13.8 Gyr, no charge/light, cold early, nearly collisionless, NOT counted in Big-Bang element-making (deuterium),
amount 5.4x ordinary matter. Limits marked (mem) are from memory of the literature."""
import numpy as np
def t_of_T(T, g): return 2.42/np.sqrt(g)/T**2
print(f"quark era ends (quarks bind into protons/neutrons, T ~ 150 MeV): t ~ {t_of_T(150, 61.75):.1e} s")
print(f"neutrinos stop trading heat (T ~ 1 MeV): t ~ {t_of_T(1, 10.75):.1f} s;  element-making (T ~ 0.07 MeV): t ~ {t_of_T(0.07, 3.36)/60:.0f} min")
sum_m = 0.064                                          # eV, DESI DR2 + CMB upper limit (mem)
print(f"neutrinos: max share of the missing mass = (sum m / 93 eV) / 0.120 = {sum_m/93/0.120:.1%}  (and they are hot)")
inv = [("photons, gravitational waves", "move at light speed: hot", "no"),
       ("quarks, gluons", "bound into protons/neutrons by ~1e-5 s -> counted by deuterium", "no"),
       ("electrons, muons, taus", "charged; muons/taus decay in microseconds", "no"),
       ("W, Z, Higgs, top quark", "decay in ~1e-25 s", "no"),
       ("neutrinos", "stable, dark, but hot and at most ~0.6% of the need", "no"),
       ("protons, neutrons, nuclei, atoms", "ordinary matter: deuterium caps it (x6 short)", "no"),
       ("pions, kaons and other hadrons", "decay in < 1e-7 s", "no"),
       ("QUARK NUGGETS (dense lumps of quarks formed at ~1e-5 s, Witten 1984)", "locked away BEFORE element-making, so deuterium never 'sees' them; cold, dark, heavy", "OPEN (needs strange-quark matter to be stable; survival in the hot early universe debated, mem)"),
       ("PRIMORDIAL BLACK HOLES", "known physics; asteroid-mass window open (mem)", "OPEN (needs big small-scale ripples)"),
       ("SPACETIME 'CLOCK DUST' (preferred-time gravity; our grid's clock)", "no particle at all; cold and collisionless automatically", "OPEN (amount not derived)")]
print("\ningredient | why | verdict")
for a, b, c in inv: print(f"{a} | {b} | {c}")
