# Ideas ledger: every idea, grouped by family, so nothing is re-proposed under a new name
Check here before suggesting anything. A new idea must say which family it belongs to and what is different from the closed tests.

## Family 1. MOND from link mechanics (links stretch, slacken, squish, change length/width, take up slack)
Closed. Uniform loose links (Saturn drift 30-40x too big), connected mosaic (snaps all at once, no square-root law), gap strengths
from geometry (slope 0.34-0.38 vs 0.56), channel geometry for a0 (needs coupling 40-900). `grid_models/`, `fluid/a0_from_channels.py`.
Includes "deeper wells stretch links more, so links respond more strongly" (EMOND-like). Do not re-propose without a new
mechanism for WHY the response changes. The math of EMOND (a0 depends on well depth) is untested, but its grid story is this family.

## Family 2. Relativistic base
AeST: closed (fails Solar System by 1e3-1e4; early instability). Khronon: passes (current base). `lorentz/`, `aest_upgrade/`.

## Family 3. Dark-energy law rho_DE ~ adot^-1/2
Passes (beats Lambda, Delta chi2 -7). Consequences: fate (`future/`), z = 0.65 release. Hubble tension: not fixed by late physics.

## Family 4. The CMB needs a cold dark component; galaxies must not end up with it; clusters need some
The core open problem. Root cause found every time: the fluid is cold and clumpy early (CMB, Lyman-alpha), galaxies form early
and capture it, and lensing shows no halos at z 0.3-0.9.
- Keep it out by pressure or stiffness: mass term, Jeans window, density pressure, DBI phase change, viscoelastic: closed.
- Keep it out by a single switch (rate, speed, stress, mass, depth): closed (galaxies in the middle).
- Add it next to MOND (hybrid, salt rule): closed.
- Move it with currents or drag (gravity currents, frame-locked, partially locked): closed or costly (growth damping).
- Other matter (hot DM, missing baryons): closed.
- Heat it out (shaken ocean, isolated wells, dark fountain): closed (lensing vs redshift; Lyman-alpha vs escape).
- The fluid in a galaxy IS its MOND halo (one thing, not added), the original ocean picture: budget checks PASS
  (supply ample; total MOND halo mass of the universe = 0.3-1.5 x the CMB fluid for g_e = 0.01-0.05 a0). Open: conversion and
  saturation mechanism; clusters still 35% short. `fluid/fluid_is_halo.*`

## Family 5. Black holes
Planck stars: consistent, no observable effect for >1e26 yr. Ended holes can't make voids or fix H0: closed.
