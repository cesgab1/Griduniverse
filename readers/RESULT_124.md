# RESULT 124: does the lock leak?
Prereg: PREREG_124.md. Code: iter124_lock_leak.py. Numbers: RESULT_124_numbers.txt. Sources: SOURCES_124_126.md.
Model: leak per interaction P = (l / lambda)^k, lambda = hbar/(mc) (reading footprint).
| Data (90% CL) | limit | k = 1: cell l < | k = 2: cell l < |
|---|---|---|---|
| electrons, VIP-2 2022 (headline) | 6.8e-42 | 2.6e-54 m | 1.0e-33 m |
| electrons, VIP-2 2024 | 6.7e-43 | 2.6e-55 m | 3.2e-34 m |
| electrons, Majorana 2023 (closed system, model-dependent) | 1.0e-48 | 3.9e-61 m | 3.9e-37 m |
| nucleons, Borexino 2010 (headline) | 4.1e-60 | 8.6e-76 m | 4.3e-46 m |
| nucleons, Borexino 2026 (read via summary page) | 7.0e-61 | 1.5e-76 m | 1.8e-46 m |
| photons forced to lock, English+2010 (wavelength ~1 um, approx.) | 4.0e-11 | 4e-17 m | 6e-12 m |
Expectations: electrons k=1 far below Planck size: HIT. k=2 electrons weak (1e-33 m, ~60x Planck): HIT.
Nucleons exclude k=2 at Planck-size cells (need < 4e-46 m): HIT.
Steepness needed for the leak to hide with real cells: nucleons require k >= 3.2 (Planck-size cells) or
k >= 5.4 (1.7e-27 m cells).

LESSON 124: the lock does not leak in any simple way. Exclusion behaves like an exact rule, not a traffic rule that
fails now and then. In the reader picture, the lock must be built into the SHAPE of a fermion's pattern (the
sign flip when two identical patterns swap, psi - psi = 0), not enforced by a busy signal on the cell. Carried to
125: capacity effects, if any, must be extremely weak -> expect no capacity correction to h/m.
