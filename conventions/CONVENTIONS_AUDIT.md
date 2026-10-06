# Conventions audit -- gravity, space-time and cosmology (Coalesce)
Rule: a convention is harmless only if varying it cannot change any prediction. Every entry is sorted into
A (pure bookkeeping), B (modelling choice that CAN change results -> must be varied), C (assumption dressed as a
convention -> is a testable hypothesis). 'Drill' = what we did or will do. (m) = from memory, to verify before use.

## A. Pure bookkeeping -- cannot change any measurement
| convention | origin | why harmless |
| units (metre, second, kg; c and h fixed by definition) | SI 1983 / 2019 | change units, predictions convert exactly |
| metric sign (-+++ vs +---), index placement | Einstein/Minkowski | notation only |
| coordinates and gauge (Newtonian vs synchronous gauge in perturbations) | GR | observables are gauge-independent; a gauge only moves bookkeeping |
| one-way speed of light (Einstein synchronisation of distant clocks) | Einstein 1905; Reichenbach | only round-trip speed is measurable; any synchronisation gives the same observations |
| Lambda written on the geometry side or the matter side | Einstein 1917 | same equations -- BUT the choice frames the 1e120 'vacuum energy' problem (drill: see C) |
| 8 pi G, h vs hbar, 'reduced' Planck mass | Planck 1899 / cosmology | factors of a few in the Planck scale; physics unchanged |

## B. Modelling choices that CAN change conclusions (vary them)
| convention | used in | effect | drill / status |
| r500 boundary (500 x critical density) | cluster masses | ties volume to total mass | DONE it.119: boundary-free points -> conclusion holds, sharper |
| mass-to-light ratio, stellar birth-mass mix (IMF: Salpeter/Kroupa/Chabrier differ x~1.6) | galaxy and cluster star masses | moves where the boost switches on | DONE drill 1: pattern survives 0.2-0.8; switch-on shifts; > ~0.6 gives impossible boost < 1 in inner points |
| hydrostatic equilibrium (gas at rest) | cluster total masses | known ~10-20% low vs lensing (m) | TODO drill: correct by lensing ratio; Coma: our 0.7-0.76 vs Tamosiunas '2x high' |
| spherical symmetry | cluster/EG masses | flattened or merging systems mis-measured | partly: Bullet done with 3-D maps (it.102/114) |
| interpolation function nu in the galaxy rule | MOND / RAR | changes cluster offsets by tens of % | TODO: rerun it.119 with 2-3 standard nu's |
| dark-energy parametrisation w0-wa (CPL, 2001) | DESI 'evolving dark energy' claim | a straight line in scale factor imposed on w | TODO drill: refit DESI+SN with other shapes (our law is one); does the 'evolution' signal survive? |
| supernova standardisation (Tripp alpha, beta; host-mass step) | dark-energy evidence, H0 | Pantheon+ vs Union3 differ | partly: hemisphere hint vanished across catalogues (it.106) |
| BAO fiducial cosmology | DESI distances | assumes a model to compress data | low risk (tested small by DESI) (m) |
| distance-ladder anchors (Cepheid vs TRGB) | H0 'tension' | ~1-2 km/s/Mpc shifts (m) | TODO if H0 becomes central |
| NFW halo shape | dark-matter fits | assumed profile shape | avoided: we use non-parametric enclosed masses |
| Bayesian priors (e.g. flat space imposed) | cosmological fits | can manufacture or hide signals | rule: report prior choice; vary when it matters |
| critical density / 'Omega' bookkeeping | cosmology | defined from H0 -- inherits H0 choices | track |

## C. Assumptions that look like conventions but are testable hypotheses
| assumption | what rests on it | test status |
| cosmological principle (same everywhere on large scales) | FLRW distances, all cosmology | tested: dipole/hemisphere/axis (it.103-106) -- no robust violation; anomaly axis p 2.6-8% open |
| averaging: smooth model = average of lumpy universe | FLRW | tested: timescape (it.99-100) disfavoured |
| G constant in space and time | Planck scale, BBN, all orbits | lunar ranging dG/dt/G < ~1e-13 per yr (m); never tested below ~50 micrometres |
| simplest Einstein action (curvature to first power only) | GR itself | Occam's choice; R^2 term fits inflation data well (Starobinsky) (m) -- open |
| matter = perfect fluid (borrowed from fluid dynamics) | cosmology T_munu | approximation; fails in collisions (Bullet) by design |
| energy conditions (energy/pressure 'reasonable') | singularity theorems -> 'Big Bang began at a point' | dark energy VIOLATES the strong energy condition -> singularity not guaranteed; supports bounce (drill: write up) |
| vacuum energy gravitates like Lambda | 1e120 problem | untested assumption; the convention in A frames it |
| redshift = expansion | all distances | tested: supernova time dilation (1+z) (m: DES ~0.5% level), CMB temperature rising with z |
| equivalence principle | GR | tested to 1e-15 (MICROSCOPE) |
| Lorentz invariance at high energy | GR, grid | tested beyond Planck energy for linear effects (GRB 090510, 221009A; our it.) |
| Gaussian, adiabatic initial ripples | CMB fits | tested: non-Gaussianity ~0 (Planck) (m) |
| Planck length = quantum-gravity scale | grid cell size | ASSUMED -- see ledger; untested |
