# The equations behind every figure (names as in GLOSSARY.md)

Every figure is backed by the equations below; the scripts that draw them (`paper/figures/make_*.py`, `orbits.py`) evaluate
exactly these. Status tags as in the glossary: NEW, BORROWED, PICTURE, OPEN.

## E1. The Grid and its tension links (Fig. 1, 2, 10)  [BORROWED form, computed on our random Mosaic]
Tension link weight between neighbouring nodes i and j (wall area over link length):

    w_ij = A_ij / d_ij

Balance of the grid around a mass (the discrete field equation; m_i = mass in cell i):

    Σ_j w_ij (φ_j − φ_i) = 4πG m_i          →  continuum limit:  ∇²φ = 4πG ρ

Draw-in (the pull, i.e. gravity):  g = −∇φ = −G M r̂ / r².
Computed on a 24,000-knot random 3-D Mosaic: pull ∝ r^(−1.96) (Newton: −2; the gap shrinks with more knots). [Fig. 2]

## E2. Tick rate and the bending of paths (Fig. 3, Fig. 0)  [BORROWED: weak-field GR, reproduced by the Khronon base]

    dτ/dt ≈ 1 + φ/c²              (clocks run slower where φ is lower: near mass)
    light: refractive index n = 1 − 2φ/c² = 1 + 2GM/(r c²);  deflection α = 4GM/(b c²)   (1.75″ for the Sun)

Ray tracing in Fig. 0 and 3 integrates  d/ds (n dx/ds) = ∇n.

## E3. Orbits (Fig. 0, 11; ν Octantis and Pluto in predictions/real_orbits.*)  [BORROWED]

    ẍ = −∇φ
    (a) lone mass, φ = −GM/r:  r = a(1−e²)/(1 + e cos θ),  T² = 4π² a³ / (GM)     (closed ellipse)
    (b) inside a Pool:  M(r) = M + M_h r³/(r² + r_c²)^(3/2),  ẍ = −G M(r) r̂ / r²  (rosette: does not close)
    (c) near two masses:  φ = −Σ_k G m_k / |x − x_k(t)|;  nearby paths separate as δ(t) ≈ δ₀ e^(λt)  (chaos)
    GR correction (not drawn): perihelion advance Δϖ = 6πGM / (c² a (1−e²)) per orbit (43″/century for Mercury)

## E4. Claim 1, the Tension-Rate Law (Fig. 4, 6, 12b)  [NEW]
Grid tension (dark energy) as a function of the Stretch rate ȧ:

    ρ_DE ∝ ȧ^(−β),  β = 1/2

Equivalent form used in all codes (q = deceleration parameter):

    d ln ρ_DE / d ln a = β q,      q = −ä a / ȧ²
    simplified (used in all fits):  q ≈ (ρ_m/2 + ρ_r − ρ_DE) / ρ_tot          (treats dark energy as w = −1 inside q)
    exact (self-consistent):        q = (ρ_m/2 + ρ_r − ρ_DE) / (ρ_tot + β ρ_DE/2)
    effect of using the exact form: Δχ² vs Λ changes by +0.1 to +0.55 (−5.4 / −6.9 / −6.3); w0 −0.91 → −0.92; the Turnover
    Point (q = 0) is identical in both forms
    w = −1 − β q / 3
    Friedmann:  H² = (8πG/3)(ρ_m + ρ_r + ρ_DE)

Consequences (β = 1/2, Ω_m = 0.31): the Turnover Point, w = −1 exactly where q = 0 (z = 0.68); w0 = −0.91; matter era
w → −1 − 1/12; ρ_DE peaks at 1.074 × today's at the Turnover Point; far future a ∝ t⁵. Two-number mimic (w0, wa) = (−0.90, −0.25).
Data: beats Λ by Δχ² = 5.4 / 7.2 / 6.8 (DESI DR2 + Planck + Pantheon+ / DES-Dovekie / Union3); free fit β = 0.63 ± 0.22.
Microscopic motivation (beta_derivation.py): independent ± exchanges across each link, N ∝ 1/ȧ per stretch time, net
tension ∝ √N ∝ ȧ^(−1/2). This is a motivation, not a derivation.

## E5. Claim 2, the Electron-Tension Transfer (Fig. 6, 9)  [NEW]

    m_e(z) = m_e,0 (1 + δ)  before the switch (z > z_s),   m_e,0 after;   δ = 0.004–0.011,  z_s ≈ 100–200
    energy handed to the grid:  ΔE = δ m_e c² per electron  →  added to the Grid tension
    recombination physics scaled: binding energies ∝ α² m_e, Thomson cross-section ∝ α²/m_e²   (camb_vc patch)
    21-cm line: ν₂₁ ∝ m_e²  →  Transfer Step at ν_s = 1420.4 MHz / (1 + z_s) = 7–14 MHz, |ΔT_b| ≤ 0.7 mK

Data: δ = 0.0099 ± 0.0049 (Planck + ACT DR6 + SPT-3G D1 + DESI DR2 + DES-Dovekie), Fig. 9.

## E6. Elasticity: stiffness, tension, squeeze and stretch limits (Fig. 10, 12)
Stiffness of a tension link [BORROWED scale; NEW reading]:

    link tension scale  F = c⁴/G = 1.21 × 10⁴⁴ N  (Planck force);   stiffness of space in Einstein's equations c⁴/(8πG) = 4.8 × 10⁴² N
    allowed drift: Ġ/G = (7.1 ± 7.6) × 10⁻¹⁴ yr⁻¹ (lunar laser ranging)  →  |ΔG/G| < 3.1 × 10⁻³ (2σ) over 13.8 Gyr;
    nucleosynthesis: |ΔG/G| ≲ 0.1 in the first minutes

Cosmic tension (energy per volume = tension per area per length):

    today  ρ_DE c² = 5.3 × 10⁻¹⁰ J/m³;   MAX 1.074 × today's at the Turnover Point;  0.22 × today's at the CMB;  → 0 in the far future

Squeeze limit (MAX compression), the Density cap [BORROWED from loop quantum cosmology, not measured]:

    ρ_cap ≈ 0.41 ρ_P,   ρ_P = c⁵/(ħ G²)   →   ρ_cap = 2.1 × 10⁹⁶ kg/m³
    collapse/expansion with the cap:  H² = (8πG/3) ρ (1 − ρ/ρ_cap)     (Bounce; Planck star; Fig. 8)

Cell size window (MIN / MAX):

    MAX  ℓ_max = √12 ħc / E_QG,2 = 5.7 × 10⁻²⁸ m,  from E_QG,2 > 1.2 × 10¹² GeV (LHAASO, GRB 221009A)   [measured]
    MIN  ℓ_P = √(ħG/c³) = 1.6 × 10⁻³⁵ m                                                                [model]

Stretch limit [OPEN requirement]:

    stretch since the Planck era = T_P / T_CMB = 5.2 × 10³¹;  allowed stretch per cell ≤ ℓ_max/ℓ_P = 3.5 × 10⁷
    → a cell that only stretched would now be 8 × 10⁻⁴ m (24 powers of ten too big): the Grid must ADD cells as space expands.

## E7. The Ocean (Fig. 5)  [BORROWED: cold dark matter]
Collisionless, pressureless fluid: ∂ρ/∂t + ∇·(ρv) = 0,  Dv/Dt = −∇φ, sourcing ∇²φ = 4πG(ρ_b + ρ_Ocean).
Pools (halos) fitted per galaxy: cored profile ρ = ρ₀ r₀³ / ((r + r₀)(r² + r₀²)) (Burkert) in Fig. 5b.

## E8. The Two-layer electron (Fig. 7)  [BORROWED construction: Kaplan domain-wall fermion]

    H = sin k σ_x + sin k_s σ_y + (M(s) − 2 + cos k + cos k_s) σ_z,   M(s) = +M inside the slab, −M outside
    light mass ∝ (1 − M)^W  (halves per layer at M = 1/2; W = layers between the two halves)

## E9. Memory version of Claim 1 (toy grid; pre-registered as the second P1 curve)  [NEW, memory DERIVED]
Jostled tension tau (velocity of a minimally coupled grid field) with Hubble friction kappa = 3 and thermal kicks:
    dX/d ln a = -2 kappa X + 1/(a H),    rho_DE = A X,   kappa = 3 (derived),  A fixed by today's dark energy
Crosses w = -1 at z = 0.46; w0 = -0.90.

## E10. Size of dark energy: the conserved stamp  [CANDIDATE, not established]
Promote Claim 1's coefficient C to a conserved quantity (unimodular-type Lagrange multiplier tau):
    S  ⊃  -∫ dt C(t) [ N a^3 (a H t_P)^(-1/2) - d tau/dt ]   ->   dC/dt = 0   (no jostling can change it)
Its value is fixed only by a cell count:  rho_DE ~ rho_P / sqrt(N4),  N4 = number of Planck 4-cells counted.
    past light cone  -> tracks the horizon: EXCLUDED (CMB; low-z leftover test, iteration 37)
    whole history    -> fixed fee, right size if N_total ~ 4.4 x the visible past: requires a finite universe (Branch A),
                        testable half = flat wrap-around space, cube side 27.5-38 Gpc (iterations 32, 38)

## E11. Primordial ripples ('watermark') from the grid's smallest-scale scaling  [BORROWED form, one FITTED number]
    omega^2 = p^(2z) / M^(2z-2),  z = 3 - eta   ->   n_s - 1 = -2 eta / (1 - eta),   eta = 0.013-0.017 (fitted; layers give
    the right order c/(16 pi^2), coefficient and sign not computed)
    general era (iteration 51): n_s - 1 = -(eta/3) eps/(1 - eps/3), alpha_s = 0; crossover route EXCLUDED (alpha_s = 8(n_s - 1))

## E12. Grid growth  [PICTURE, tested]
Links stretch with expansion (strain += d ln a), break at a tolerance, a new cell is inserted: explains cell addition and the
sqrt(N) fluctuations; does NOT set the memory or the size (iteration 34).

## E13. Gravity on the random Mosaic: curvature read by neighbourhood averaging  [BORROWED definition; tested on our grid]
On each time slice (preferred slicing), curvature at grid point x is read from how its neighbourhood moves onto a neighbour's:
    kappa_eps(x, y) = 1 - W(m_x, m_y) / d(x, y)          (Ollivier-Ricci; m_x = uniform over the points within distance eps of x,
                                                           W = cheapest transport cost between the two neighbourhoods,
                                                           costs = link lengths)
    R_eps(x) = (2 n (n + 2) / eps^2) < kappa_eps(x, y) >_y   ->   n = 3:  R_eps(x) = (30 / eps^2) < kappa_eps(x, y) >_y
Grid gravity action (ADM form, preferred slicing; Claim 1 term as in E4):
    S = (1/16 pi G) SUM_t dt SUM_x V_x N [ K_ij K^ij - lambda K^2 + R_eps(x) ]  +  S_Claim1  +  S_matter
    V_x = proper volume per point (1/density for a sprinkling uniform in proper volume)
Scale window:  cell spacing a  <<  eps  <<  curvature radius;  at least ~1000-2000 points per neighbourhood (eps >~ 6-8 a).
WHY (iterations 39, 41, 41b):
  - A random grid has 'sliver' cells. Regge's definition reads curvature from the angles around ONE link, so a sliver gives a
    wildly wrong angle and these errors never average away: the Regge action on random grids did not converge (39).
  - Reading curvature from how whole NEIGHBOURHOODS spread or crowd averages over thousands of points, so the randomness
    cancels: on the random grid it gives the true curvature (1.18 +/- 0.12 at ~2100 points per neighbourhood, 41b).
  - Same principle as the rest of the model: the layers average the jostling; the geometry is read by averaging too.
  - Consequence: gravity has a built-in resolution eps (several cell spacings, i.e. ~1e-31 to 1e-27 m): below eps the grid has
    no well-defined curvature. Far below any experiment; it is the grid's own 'pixel size' for gravity.
Status: spatial curvature term TESTED on the random grid; the time-derivative part (K_ij) on a random grid NOT yet tested; the
full Einstein-action convergence from R_eps NOT yet shown (only for tidy cells with Regge, E12/iteration 39b).

## E14. Two levels of the grid (Coalesce's framing: tidy for GR, different at the quantum level)
    Level 1 (quantum): random Mosaic cells, spacing a = 0.49 sqrt(N) l_P (E/iteration 40); no curvature per cell.
    Level 2 (classical): neighbourhoods of radius eps >~ 6-8 a (~1000-2000 cells), curvature R_eps (E13) -> Einstein gravity.
    Switch-over scale eps: below it the grid looks different (random; effective dimension expected to fall towards ~2, as in CDT
    and in Horava scaling with z ~ 3: d_s = 1 + 3/z); above it, smooth and 'tidy'.
WHY: general relativity only ever sees averaged neighbourhoods, never single quantum cells. Tested: tidy cells give Einstein's
action (39b); averaging random cells gives the right curvature (41b). Not yet shown: the full action built from averaged random
cells, and the effective dimension of our own grid below eps.

## E15. The averaged action IS Einstein's (spatial part; pen-and-paper, quantum_gravity/emergence/averaged_action_derivation.md)
    S_eps = (1/16 pi G) INT sqrt(g) [ R + eps^2 (a R^2 + b Ric_ij Ric^ij) + O(eps^3) ] + boundary terms
From: Ollivier's expansion (borrowed) + six-direction averaging (trace -> R, +/- pairs cancel gradient terms) + total derivatives
integrate away. Leading coefficient measured on our grid (0.99 +/- 0.08 at the bump centre). Corrections in the real universe:
curvature-squared <~ 2e-61 relative (L ~ 10 km), random noise ~ 3e-39 or smaller per cm^3. a, b not computed; time part not yet.

## E16. Time part on the random grid (iteration 45)  [TESTED]
    K_ij fitted per point from link-length changes between ticks:  d ln(l_e)/dt = n_e . K . n_e   (least squares over ~34 links)
    kinetic density K_ij K^ij - lambda K^2;  uniform expansion -> -6 H^2 (exact); gravitational wave -> A^2 w^2/4 (0.993)
Together with E15 (spatial part): the averaged random grid gives the full ADM action
    S = (1/16 pi G) INT dt d^3x N sqrt(g) [ K_ij K^ij - lambda K^2 + R + eps^2 (curvature^2) ]  (+ Claim 1 term, matter)
which is Einstein's gravity for lambda = 1 (Horava-type for lambda != 1; lambda is a parameter of the base, not derived here).

## E17. Two sectors (quantum_gravity/lambda/two_sector.md)  [STRUCTURE; leakage factor c not computed]
    Gravity: N layers on a grid random in space AND time -> frame-free -> lambda = 1; small-scale dimension -> 2 (causal-set result).
    One frame-dependent field (tension/slicing): carries Claim 1, the 'now' (constant-K slices) and the z = 3 watermark.
    Leakage into gravity ~ c/N.  Graviton speed (GW170817, 1e-15) -> N >~ c x 1e15; for c = 1: N ~ 1e15 - 5e15,
    cell 2.5 - 5.7 e-28 m (at today's gamma-ray-burst limit).

## E18. The equations run on the grid (iterations 47-49)  [MEASURED / DERIVED / BORROWED as marked]
    history: N_n found from dS_grid/dN_n = 0 -> H(z) = Friedmann + Claim 1 to 2e-9; acceleration equation holds unimposed (~step^2)  [MEASURED]
    ripples: d2h_i/dt2 = (1/norm) SUM_j exp(-r_ij^2/2s^2)(h_j - h_i) on random points -> speed = continuum to 0.03%, isotropic;
             coherent loss per distance ~ k^2.7 (real universe: ~1e-40 over 40 Mpc)                                            [MEASURED]
    leakage: tree level 0 (VCDM class) ; Delta lambda = Omega_DE/6 on horizon scales only ; K in Claim 1 averaged over > 0.03 mm [BORROWED/DERIVED]
    electron switch: m_e = 1 + delta/2 [1 + tanh((z - z_s)/dz)] requires z_s/dz >~ 5 (methanol, z = 0.89)                     [MEASURED bound]

## E19. Growth of structure with the frame field (iteration 50)  [BORROWED form (VCDM), MEASURED numbers]
    Phi' + aH Phi = (3 a^2/2k^2) rho_m theta F(k),   F = [k^2 - 3 a^2 dH/dt]/[k^2 + (9/2) a^2 rho_m]
    our model vs same history with ordinary gravity: < 1e-5 for k >= 0.01 h/Mpc; 1.2% at k = 3e-4; ISW -3.6% at k = 2e-4
    -> unobservable (cosmic variance); background-only fits justified.

## E15 status note (iteration 52)
The correction term eps^2 (a R^2 + b Ric^2) is an ASSUMED form: a direct measurement of a, b failed (numerics not converged),
and Ollivier's theorem allows an eps^3 term in kappa (i.e. eps^1 in R_eps). Leading term (Einstein) confirmed to 0.6%.
Any form is unmeasurable (<= 5e-31 relative at the most curved tested places).

## E20. Size of dark energy from the QG piece (iteration 53)  [BORROWED core (Sorkin), DERIVED with our cells; 3 free choices]
    rho_DE = c_rule sqrt(N) hbar c / (a^2 sqrt(V_book)),   a = 0.49 sqrt(N) l_P   ->   rho_DE = c_rule (4.15/sqrt(N)) rho_P / sqrt(N4_book)
    whole book: S x T = 4.0 c_rule^2 x 17.2 / N ;  space >= CMB topology floor (S >= 1.70)  ->  N <= 41 c_rule^2
    consequences: cells 1-6 l_P, cutoff ~1e19 GeV (no photon delays); finite flat wrap-around space just beyond 27.5 Gpc; time ends

## E21. Size of dark energy, non-quantum route: the electron pays (iteration 54)  [NEW; 2 free choices: criterion, one-way]
    rho_DE(z_s) = delta m_e c^2 n_e(z_s),  z_s: Compton heating rate = expansion rate (z_s = 124, CAMB)
    -> delta predicted 0.60% (Claim 1) vs measured 0.99 +/- 0.49%;  handover must be one-way (reionisation re-couples the gas)
    fork with E20: electron route (delta = 0.6%, 21-cm step 11.4 MHz) vs cell-tally route (finite wrap-around space)

## E22. Two-level picture of dark energy's size (iteration 55)  [BORROWED structures: unimodular non-conservation, sequestering]
    QG (how much):  rho_DE = c_rule (4.15/sqrt N) rho_P / sqrt(N4_book)
    GR (who pays):  d(stamp)/dt = J >= 0 (one-way, a tally),  J = electron mass-energy released at z_s = 124  ->  delta = 0.60%
    GR twin of the tally: (1/4)<rho_m>_book = 0.18-0.41 rho_DE for books ending 0-5 Gyr from now (automatic, 'why now')

## E23. Electron-pays mechanism (iterations 56-60)  [LEAD; exact identity + motivated rule]
    open pop: e^2 / (4 pi eps0 (hbar/m_e c)) = alpha m_e c^2   ->   delta = alpha = 0.730%
    only point-like charges hold pops (size < hbar/mc): electrons yes; protons, nuclei no
    rho_DE(z_s) = alpha m_e c^2 n_e(z_s),  z_s ~ 116 (end of the gas's thermal tie to the light), then Claim 1
