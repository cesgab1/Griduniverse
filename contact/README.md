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
Status: candidate, not tested on galaxies. Link to existing work: the electron switch (z ~ 100-200) is a moment where the
grid already exchanges energy with matter; a contact that turns on then would avoid the CMB bound but hit 21-cm.
