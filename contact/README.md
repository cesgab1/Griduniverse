# Family 7: the fluid in contact with visible matter (Oct 2026)
Why: without MOND the one open issue is that galaxies are tighter and more regular than halos placed by abundance matching
(no_mond R1-R3). A halo that trades energy with the visible matter it surrounds would settle to match it.
Literature version: Famaey, Khoury & Penco 2018 (arXiv:1712.01316) show DM-baryon energy exchange gives the radial
acceleration relation if sigma = C a0/(n_DM v^2), C ~ 1/16; clusters never equilibrate (so no MDAR there, as observed);
BTFR normalisation needs simulations. Follow-up (Baryon-Interacting DM, OSTI 1802503) claims BTFR + central surface density.
Grid reading: each bit of visible matter trades energy with the fluid at a fixed rate per baryon, set by a0, however much
fluid is there. That is an environment-dependent coupling; a0 is still put in by hand (a0 ~ cH0/6 is the only hint).
Checks (contact_bounds.*):
- needed in galaxies: sigma/m ~ 0.3-18 cm^2/g;
- the same law at recombination: 6-110x above the Planck bound on DM-baryon scattering (rough: the bound is for constant
  sigma) -> the contact must switch on after recombination;
- intergalactic gas: exchange faster than expansion at z = 3 (Gamma/H ~ 7), ~0.6 at z = 100 -> must be checked against the
  Lyman-alpha forest temperature and the 21-cm signal (z ~ 15-20).
Follow-up (early_tests.*), heavy dark matter as the mechanism needs:
- 21-cm: PASSES for m >= 100 GeV (the cold fluid can cool the gas by <= 5%); light (1-10 GeV) would deepen the signal 33-83%.
- Lyman-alpha: FAILS. The same law drags the intergalactic gas onto the fluid in ~50 Myr, while pressure needs 1.2-2.2 Gyr to
  smooth it over the measured filtering scale (~80 kpc): the gas would be locked 25-43x too fast and the forest would show no
  pressure smoothing (it does). Generous choice (v = gas sound speed); the real fluid spread is smaller, making it worse.
- CMB: 6-110x over Planck (above), unless switched on after recombination.
- Inside galaxies: a scattering contact also drags. Disc gas rotating through a non-rotating halo spins down in 0.3-0.55 Gyr
  (discs are ~10 Gyr old; gas and stars co-rotate). Routing the momentum into the grid instead makes it a drag toward the grid's
  rest frame, which is worse.
Verdict: FAILS. Closed as a universal law (Lyman-alpha, CMB) and as a galaxy-only law (disc gas spin-down). The open R2/R3
issue stays where standard LCDM has it: galaxy formation physics.
Old status line: candidate, not tested on galaxies. Link to existing work: the electron switch (z ~ 100-200) is a moment where the
grid already exchanges energy with matter; a contact that turns on then would avoid the CMB bound but hit 21-cm.
