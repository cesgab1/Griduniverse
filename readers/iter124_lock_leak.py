"""PREREG 124: lock-leak model P = (l / lambda)^k  ->  l < lambda * P_limit^(1/k)."""
lam_e, lam_N = 3.8616e-13, 2.1031e-16          # reduced Compton wavelengths hbar/(mc), m (electron, nucleon)
lam_ph = 1e-6                                  # photon wavelength in Bose-symmetry test, ~1 micron (approximate)
tests = [  # (label, footprint, limit, headline?)
 ("electron, VIP-2 2022 open system (90%)",              lam_e, 6.8e-42, True),
 ("electron, VIP-2 2024 modulated current (90%)",        lam_e, 6.74e-43, False),
 ("electron, Majorana 2023 closed system (90%, model)",  lam_e, 1.0e-48, False),
 ("nucleon, Borexino 2010 strong (90%)",                 lam_N, 4.1e-60, True),
 ("nucleon, Borexino 2026 full data (90%, via summary)", lam_N, 7.0e-61, False),
 ("photon (boson forced to lock), English+2010 (90%)",   lam_ph, 4.0e-11, True),
]
lP, lgrid = 1.616e-35, 1.7e-27
out = []; P = lambda s="": (print(s), out.append(s))
P(f"{'test':52s} {'k=1: l <':>11s} {'k=2: l <':>11s}")
for lab, lam, lim, head in tests:
    P(f"{lab:52s} {lam*lim:11.1e} {lam*lim**0.5:11.1e}" + ("  <- headline" if head else ""))
P(f"\nreference: regular-grid limit from photon timing 1.7e-27 m; Planck length {lP:.1e} m (comparison only)")
for lab, lam in (("electron", lam_e), ("nucleon", lam_N)):
    P(f"predicted leak if l = Planck length: {lab:8s} k=1 {lP/lam:.1e}, k=2 {(lP/lam)**2:.1e}")
open("RESULT_124_numbers.txt","w").write("\n".join(out)+"\n")
