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

# Same-fashion comparison of four geometries (geometries.py, compare_geometries.*, calibrate_geometries.*, systematic_checks.*)
All four in the same ball (radius 16, ~17,000 nodes), same point mass at the centre, same boundary, same solver, same scores,
two link rules (static links; links that thicken with flux, sqrt rule, saturating). Mosaic and fungal: two random seeds each.
Geometries: cubic = symmetric lattice; mosaic = random Delaunay tiling (degree 15); web = 3-D spider web, spokes + rings
around the centre; fungal = grown hyphae that branch and fuse (degree 3.2).

| geometry | static: outer slope / TF | adaptive: outer slope / TF | Newton inside (adaptive/static) | MOND shape rms | cubic-axis pattern | off-centre mass |
|---|---|---|---|---|---|---|
| cubic (symmetric) | -2.05 / 2.00 | -1.08 / 1.00 | 1.00 | 0.10 | +1.95 (strong) | -0.84 |
| mosaic | -1.92, -2.10 / 2.00 | -0.95, -1.05 / 1.00 | 1.00 | 0.05 | -0.05 (none) | -0.87 |
| web | +0.30 / - | +0.30 / - | (never Newtonian) | 0.68 | 0 (built around one centre) | -0.59 |
| fungal | -1.68, -1.78 / 2.00 | -0.83, -0.91 / 1.00-1.02 | 1.00 | 0.08 | -0.00 (none) | -0.83 |
MOND targets: outer slope -1, Tully-Fisher (TF) 1, Newton inside 1.00, rms 0. Newton: outer -2, TF 2.
(Outer slopes of ~ -0.85 for off-centre masses partly reflect the nearby boundary; all non-web cases agree there.)

Verdict:
1. With static links every geometry gives plain Newton: no geometry alone makes MOND. The link rule does it.
2. With the flow-thickening rule, mosaic, fungal and cubic all give MOND (1/r far out, Newton inside, Tully-Fisher slope 1).
3. The symmetric cubic lattice imprints its axes: a direction pattern worth ~2 grid spacings at r = 10 (~20%). The
   nonlinear rule makes the lattice's directions matter at every scale. That conflicts with round lensing. EXCLUDED.
4. The spider web is never Newtonian (fixed spokes give a pull that doesn't fall off), and it has a built-in centre, so it
   can't serve many masses. EXCLUDED.
5. Mosaic and fungal tie. Both are statistically isotropic with no pattern. Mosaic fits the MOND shape slightly better (rms
   0.05 vs 0.08) and fungal is cleaner close in. The geometry only needs to be a random, connected 3-D mesh.
   The fungal picture adds a reason for the rule itself (growth follows flow).
Earlier mosaic failures were failures of the slack rule, not of the mosaic geometry.
