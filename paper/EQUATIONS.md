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

    d ln ρ_DE / d ln a = β q,      q = −ä a / ȧ² = (ρ_m/2 + ρ_r − ρ_DE) / ρ_tot  (to first order)
    w = −1 − β q / 3
    Friedmann:  H² = (8πG/3)(ρ_m + ρ_r + ρ_DE)

Consequences (β = 1/2, Ω_m = 0.31): the Turnover Point, w = −1 exactly where q = 0 (z = 0.68); w0 = −0.91; matter era
w → −1 − 1/12; ρ_DE peaks at 1.074 × today's at the Turnover Point; far future a ∝ t⁵. Two-number mimic (w0, wa) = (−0.90, −0.25).
Data: beats Λ by Δχ² = 5.5 / 7.3 / 6.8 (DESI DR2 + Planck + Pantheon+ / DES-Dovekie / Union3); free fit β = 0.63 ± 0.22.
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
