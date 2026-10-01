# Fluid-tube network that thickens with flow (fungus / slime-mould type grid) -- adaptive_network.*
Coalesce's two ideas combined: the planks are tubes of fluid with dead space between (like fungal hyphae), and the network
grows like a fungus/slime mould: tubes carrying more flux thicken, idle ones thin (Tero et al. Physarum rule, steady D = F(|Q|)).

Why it matters: in steady state the flux obeys div(D(|grad p|) grad p) = source, the AQUAL (MOND) field equation.
- Tubes thicken as sqrt(flux) while thin, and saturate when fully built: deep MOND far out, Newton close in. a0 = the flux per
  unit area at which tubes are fully built.
- The switch radius then sits automatically where the pull falls to a0, sqrt(GM/a0). That is the condition the spider web
  needed by hand, and it gives v^4 ~ M (baryonic Tully-Fisher) for free.
- Thick near the centre, thinning outward: the "tighter near the centre" picture.
- Exponent < 1 keeps loops (a mesh, not a tree), so the network stays roughly the same in all directions and doesn't
  avalanche like the mosaic.

3-D lattice test (41^3, point source, adapt to steady state):
- Field falls as 1/r beyond the MOND radius and tracks the MOND prediction to ~5-20%. The lattice shifts the effective a0 by a
  constant factor of ~1.2.
- Along the axes vs the diagonals the field agrees to 5-10% (cubic-lattice artefact, larger near the boundary).
- Tully-Fisher: outer v^4 vs M slope 1.01 with fully-built saturation (MOND: 1). The smooth rule gives 1.19, because these
  small lattices haven't reached deep MOND.
Still put in: the sqrt reinforcement exponent (1/2) and the saturation flux a0. A linear rule (exponent 1) gives a constant pull
(no 1/r), so the exponent matters.

Open checks:
(1) Bookkeeping with Khronon (one conserved fluid = MOND halo): does the fluid held in the thickened tubes equal the MOND phantom
    mass? That needs tube conductance ~ (fluid content)^(1/2) in deep MOND, and little extra fluid inside the Newtonian core.
(2) Growth takes time: MOND would build up as a galaxy forms. z ~ 2 disks with falling rotation curves (Genzel+2017) fit
    'not yet grown'. Lensing is unchanged between z 0.3 and 0.9, so growth must take <~ 2 Gyr.
(3) Solar System: the saturation must be sharp (Cassini), as for every MOND theory.
