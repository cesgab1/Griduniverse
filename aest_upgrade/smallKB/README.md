# Can AeST's cosmology survive the small K_B the Solar System requires?

Solar System / pulsar preferred-frame bound (lorentz/aest_alpha1.*): alpha1 = -4[K_B + (2-K_B)/(1+lambda_s)] -> K_B < 2.5e-5.

Cosmology (AeST eqs 7-12 of Skordis & Zlosnik 2021, re-checked verbatim against the paper; our CLASS patch):
- fast_modes.py: the vector-sector block (alpha, E) has a growing mode in the early universe with rate ~ sqrt(3 Omega_fld/K_B) x aH,
  from the K'chi term of eq. 12 (K'Q = 8 pi G rho > 0). Independent of the shape of K(Q). It switches to harmless oscillation only once
  (H+Q)Q exceeds ~ rho, i.e. when H drops below ~Q0, which the MOND-scale requirement (mu = sqrt(2K2/(2-K_B)) Q0 <~ 1/Mpc) keeps late.
- growth_exponent.py: total amplification e^N:  see growth_exponent.txt  (K_B = 0.5: modest; K_B = 2.5e-5: astronomically large)
- sigma8_vs_KB.py (CLASS, high precision): the mode leaks into the matter spectrum through the pressure term:
  K_B 0.5 -> 0.2: P(k=1) +8%;  0.15: x2.05;  0.1: x633, sigma8 0.74 -> 1.50  (Cosh parameters).
  The paper's Exp model uses K_B = 0.1 with Z0 = 1e-17, which suppresses the leak (c_ad^2 ~ Z0/Q), but cannot stop the field itself growing.

Verdict: in its published form AeST cannot satisfy both the preferred-frame bound (K_B < 2.5e-5) and early-universe stability
(K_B >~ 0.1-0.2). Caveats: linear analysis of our implementation of the published equations; alpha2 not computed; a modified vector
sector (e.g. extra aether terms c2, c4 chosen to cancel alpha1 with c14 = 0) is the obvious repair to explore.
