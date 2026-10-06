# PREREG 113 -- emergent gravity (Verlinde 2016) in Coma, with the four kinds of entropy (Coalesce)
Committed BEFORE running iter113_emergent.py.
Idea: gravity is not fundamental; it emerges from quantum information (entanglement) of space. In de Sitter space (our
universe with dark energy) the entanglement entropy has an extra VOLUME-law part; matter displaces it, and where the
displaced entropy is small (weak pull, g < ~cH0/6) the response looks like extra 'apparent dark matter'.
Formula (spherical, static only): M_D(r)^2 = (c H0 r^2 / 6G) d[M_B(r) r]/dr. No free parameter (H0 = 70).
Data: Coma gas beta model, Briel+1992 as used by Bonafede+2010 (arXiv:1002.0594, verified): beta 0.75, r_c 291 kpc,
n_e0 3.44e-3 cm^-3 (H0 71); mu_e 1.17. Stars 5e13 Msun within 2.5 Mpc (project value), same shape as gas (sensitivity:
Hernquist a = 0.15 Mpc). 'Measured' mass: hydrostatic from the same gas, isothermal kT = 8.2 keV:
M(r) = 3 beta kT r^3 / (G mu m_p (r^2 + r_c^2)), mu 0.6; check against the project's 1.6e15 at ~2.5 Mpc.
Report M_EG/M_measured at r = 0.1, 0.3, 0.5, 1.0, 1.5, 2.5 Mpc.
Entropy ledger for Coma (R = 2.5 Mpc), each in units of k_B:
 thermodynamic (gas Sackur-Tetrode, stars, CMB, neutrinos); coarse-grained gravitational (black holes, Bekenstein-Hawking);
 entanglement (area law A/4l_P^2 = holographic maximum; de Sitter volume-law part S_DE = (A/4l_P^2)(r/L), L = c/H0);
 entropy displaced by matter S_M = 2 pi M c r / hbar; von Neumann entropy of the whole (pure state, unitary).
 Ratio S_M/S_DE = g_N/(cH0/2) is where emergent gravity switches on.
Free choices: star profile (1), isothermal T (1), H0 (fixed 70). Look-elsewhere: 6 radii, 1 model.
Expectations (memory of published cluster tests, Ettori+2017, Halenka & Miller 2018):
- outer Coma (1.5-2.5 Mpc): EG within ~30% of measured mass.
- inner (<= 0.3 Mpc): EG short by >~ x2.
- S_M/S_DE < 1 everywhere outside ~0.1 Mpc (cluster is in the 'weak pull' regime).
- Bullet Cluster: not computable (formula needs spherical symmetry); EG acts like MOND with a0 = cH0/6, so expected to
  share the QUMOND miss of iteration 102 (5 sigma). Stated, not run.
- Black-hole entropy dominates the thermodynamic/coarse-grained totals by >1e20; all far below the area-law maximum.
