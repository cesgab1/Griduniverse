# Grid Universe

A tension-grid picture of gravity, galaxies, cosmology and the electron, tested piece by piece against public data.

Space is modelled as a random grid under tension:
- mass pulls the grid into dimples, which is gravity;
- the grid goes slack where the pull is weak, giving extra pull in galaxy outskirts;
- the grid's cosmic tension is the dark energy;
- particles are knots or patterns in the grid (the electron as a two-layer pattern).

Every idea was turned into a calculation and compared with real measurements. Failures are recorded alongside successes.

`page/index.html` is the full scoreboard: what was tested, the data used, the result, and the verdict for each idea.

## Headline results

| Idea | Data | Result |
|---|---|---|
| Dark-energy law `d ln ρ_DE / d ln a = β q` (ρ_DE ∝ ȧ^(−β)) | Planck + DESI DR2 + DES/Pantheon+/Union3 | Better than ΛCDM (Δχ² −5 to −10); free fit gives β = 0.63 ± 0.22 |
| Galaxy spin from light alone (slack rule) | SPARC, 171 galaxies, blind half-split | 78% within 20% (Newton alone: 6%) |
| Disc-galaxy lensing from light alone | KiDS-1000 (Brouwer+2021) | measured/predicted 0.91–0.94 ± 0.06 |
| "Mass beats tension" boundary for slack (no free number) | KiDS-1000, 4 mass bins | fits as well as a free fade (χ² 145.5 vs 150.7 without) |
| Heavier electron at recombination (grid tension sets the constants) | Planck + ACT DR6 + DESI + DES | m_e = 1.0078 ± 0.0044, Δχ² −2.0 at +1%, H₀ ≈ 69–71; first Hubble route to survive ACT |
| Electron relaxes after recombination, its lost energy becoming the grid tension (dark energy) | Planck + ACT DR6 + DESI + DES; BBN, quasars, FIRAS, τ | Switch at recombination excluded (Δχ² +117); switch after z ≈ 900 keeps the full effect; energy matches today's dark energy if the switch is at z ≈ 100–200. Prediction: 21-cm step at 7–15 MHz |
| Forward from the Big Bang: CMB-only start, β = ½ and electron shift from the energy budget (nothing tuned to today) | Planck + ACT DR6 → predict DESI DR2 BAO, DES-Dovekie | Predicts today's distances better than ΛCDM: Δχ² ≈ −8 (BAO −3.5, SN −4.5); H₀ 67.7–69.1 |
| Upgraded AeST (AeST dark sector + ρ_DE ∝ ȧ^−½) in one Boltzmann code (our CLASS build) | Planck + ACT DR6 + DESI DR2 + DES-Dovekie | vs AeST+Λ: Δχ² −7.0 (no exchange, S₈ 0.809), −3.4 (with exchange, S₈ 0.776) |
| Early stored stretch energy | Planck + ACT DR6 | rejected by ACT (Δχ² +13 at δ = 0.03) |
| Five late-universe Hubble fixes | Pantheon+SH0ES, DESI, BBN | all ruled out (no-CMB H₀ = 68.6 ± 0.6) |

## Status (hostile review)

This is a research program, not a unified theory. Since October 2026 its relativistic home is Khronon theory (Blanchet & Skordis 2024): the grid's layers are Khronon's preferred slicing, the fluid between them is its dust-like condensate, and the MOND response sits on top. AeST, the previous base, was dropped after we found it fails the Solar System preferred-frame test and is unstable where it would pass (`lorentz/`, `aest_upgrade/smallKB/`). A dark-matter-like ingredient is still needed in the CMB fits (now supplied by the condensate) and in clusters (75–93%), but not in galaxies.

What is original and testable here:
- the dark-energy law ρ_DE ∝ ȧ^(−½);
- the forward-from-the-Big-Bang prediction (about 2σ);
- the electron-energy/dark-energy link, which predicts a 21-cm step at 5–26 MHz and a 0.4–1.1% early electron shift measurable by CMB-S4;
- several empirical findings on a₀ and clusters.

The rest is borrowed physics (credited below), tested-and-failed ideas (kept on record), or metaphor.

## Lessons learned

Every test so far, grouped by verdict. "Ours" marks results that appear not to be in the literature.

### Passed
| Test | Result | Where |
|---|---|---|
| Dark-energy law ρ_DE ∝ ȧ^(−½) on the Khronon base | beats Λ by Δχ² −7.1 (no exchange) / −3.6 (exchange), no extra parameter; −6.7 / −3.8 with a Lyman-α-safe fluid | `aest_upgrade/khronon/` |
| Forward from the Big Bang (CMB-only start) | predicts today's BAO + supernovae better than ΛCDM, Δχ² ≈ −8 | `cosmology_fits/forward/` |
| Galaxy spin from light alone (MOND rule) | 78% of SPARC galaxies within 20% | `galaxies_lensing/` |
| Disc-galaxy lensing from light alone | KiDS ratio 0.91–0.94 ± 0.06 | `galaxies_lensing/` |
| a₀ = cH₀/6 | within 1% for bulge-free SPARC at stellar M/L 0.6 (conditional on M/L) | `grid_models/taut_refill*` |
| Fluid in narrow channels → pressure ∝ density³ (the MOND superfluid law) | exact 1D Bose gas: exponent within 5% of 3 for γ > 38 | `fluid/channel_eos.*` |
| Grid seen by light | gamma-ray-burst limit needs spacing < 6e-28 m; a Planck grid passes by 7.5 orders | `lorentz/` |
| Khronon base | GR in the Solar System at 1PN, gravitational waves at c, lensing = dynamics, stable; matches ΛCDM to 0.05% with Λ | `aest_upgrade/khronon/` |
| Our CLASS AeST/Khronon build | agrees with independent CAMB run (2440.6 vs 2441.2) | `aest_upgrade/` |

### Plausible / open
| Item | Status |
|---|---|
| Energy exchange between dark energy and the fluid | no-exchange fits better on Khronon; theory doesn't yet say which |
| Electron energy → dark energy (switch at z ≈ 100–200) | consistent; predicts a 21-cm step at 7–15 MHz |
| Cassini quadrupole vs rotation curves | sharp switch passes Cassini but fits galaxies worse; shared by every MOND theory (literature 8.7σ, 1.9σ without bulges) |
| Origin of a₀ | Verlinde-type argument gives cH₀/6; not derived in the fluid picture |
| Black holes as Planck stars (no singularity; bounce/leak) | consistent with the grid; no observable consequence for astrophysical holes for >1e26 yr |
| Supermassive black-hole seeds from fluid collapse | hypothesis only (`black_holes/`) |
| **'Ocean' hybrid: MOND from baryons through the geometry + unboosted fluid halo of the cosmic ratio (5.4 × baryons, ~100 kpc core)** | one rule, R not fitted: SPARC Δχ² +26 / 2788 pts; KiDS lensing χ² 202 → 129; clusters 0.79–1.13 of measured (MOND alone 0.33–0.35). Needs a theory where only baryons source MOND and light shares their metric; core size and galaxy retention not derived (`fluid/hybrid_ocean.md`) |
| Khronon paper's DBI setting (λ_D ≈ 1) | not ΛCDM-like in our solver (σ₈ 0.22); needs λ_D ≈ 10⁷–10¹²; to be cross-checked with the authors' code if available |
| Hubble tension | treated as a calibration issue (working assumption); FRB H₀ undecided |

### Failed (and kept on record)
| Test | Why it failed | Lesson |
|---|---|---|
| Five late-universe Hubble fixes; early stored stretch energy | data (ACT Δχ² +13; no-CMB H₀ = 68.6 ± 0.6) | the Hubble tension isn't fixed by late physics in this model |
| Electron switch exactly at recombination | Δχ² +117 | timing matters; later switch survives |
| Light extra fermions for clusters | abundance 20–30× too high | cluster mass needs another route |
| Uniform loose links (slack offset) | leftover constant pull: Saturn drift 30–40 mas/yr vs ~1 allowed | check the Solar System early for any MOND mechanism |
| Connected random mosaic taking up slack | snaps taut all at once (rigidity avalanche), no square-root law | collective network behaviour differs from independent strands |
| Gap strengths from mosaic geometry | weak-pull slope 0.34–0.38 vs 0.56 measured | geometry alone doesn't give the square-root law |
| 10-µm fluid channels as the grid light travels on | excluded by 22 orders | the fluid must be dark; light rides a Planck-fine grid |
| Standalone superfluid dark matter | fails weak lensing (light ignores the phonon push; literature) | the push must act through the geometry light follows |
| a₀ from channel geometry alone | needs coupling α ≈ 40–900 | a₀ is not yet derived from the fluid |
| AeST preferred frame (ours) | α₁ = −4[K_B + (2−K_B)/(1+λ_s)]: fails by 10³–10⁴ at K_B = 0.1–0.5 | always run the PPN/preferred-frame check on a base theory before building on it |
| AeST at small K_B (ours) | early-universe vector mode grows as √(3Ω/K_B)·H: e¹³⁰⁰ at the K_B the Solar System needs | a fix in one regime can break another; check stability |
| AeST + c₄ J² repair | cancels α₁ but replaces MOND with Newton×(2/|c₄|) at low accelerations | MOND, preferred frame and stability are tied together |
| Ended black holes (stellar clusters or supermassive) creating voids/cold regions that explain the Hubble tension | none has ended (shortest lifetime >1e26 yr); all black holes together hold ~6e-4 of the energy needed; cold voids slow visible light by ~1e-29, and H₀ uses brightness and redshift, not travel time | `black_holes/ended_holes_hubble.*`; the local-void route was already excluded by supernovae |
| Khronon mass term (μ) supplying cluster mass | enough mass for cluster cores (1/μ ≈ 1 Mpc) multiplies isolated-galaxy lensing at 0.3–1 Mpc by 5–500×; for 1/μ ≲ 0.5 Mpc no static solution exists | `clusters/khronon_mass_term.*` |
| Dark fluid clumping only in clusters (Jeans/crossover length 1–4 Mpc) | removes 15–60% of small-scale power at z ≈ 3; Lyman-α forest allows ~2% | `clusters/condensate_jeans.*` (also: the AeST paper's fluid setting used in our fits is Lyman-α-excluded; with a safe setting the dark-energy result still holds, −6.7 / −3.8) |
| Fluid pressure (power law in density) keeping it out of galaxy halos but not clusters | capped by the early universe: largest Jeans mass at halo density 3e5 Msun (need >1e12). A saturating DBI pressure escapes this (see Plausible) | `clusters/jeans_window.*` |
| DBI phase change: galaxies halo-free, clusters dark-matter-like | works today (static and time-dependent), but galaxies form when the threshold was 10–90 km/s, capture halos, and those halos are self-sustaining; at KiDS lens redshifts MOND + retained halo gives χ² 50,514 vs 202 for MOND alone | `clusters/README.md` §7–10. General lesson: a cold cosmological dark fluid in AeST/Khronon-type theories must not build galaxy halos, and nothing found so far prevents it |
| Rate-dependent stiffness (viscoelastic fluid) keeping it out of galaxies | galaxy and cluster collapse times differ only 1.5× (0.84–0.94 vs 1.17–1.30 Gyr); best case still gives halos to 15–25% of galaxies | `fluid/state_change_tests.md` |
| Any single switch (rate, speed, stress, mass, potential depth) | galaxies sit in the middle of every axis: some structures that must clump lie below them and some above | `fluid/trigger_axes.*` |
| Raw "everpresent Λ" random walk | BAO + SN Δχ² ≥ +1100 | the √N idea works smoothed (β = ½), not literally |
| Wang–Unruh hidden vacuum energy | needs a cutoff tuned inside an exponential; predicts w = −1 | moving fine-tuning around is not a derivation |

### Overall lessons
1. Check the hard tests first: Solar System, gravitational-wave speed, lensing, stability. We spent weeks on AeST before finding it fails the Solar System.
2. Collective behaviour beats intuition about single parts: a single loose link isn't a network of loose links.
3. Derived numbers are often conditional: a₀ = cH₀/6 holds only for heavier stellar mass-to-light ratios.
4. The literature often already has the failure (superfluid lensing); search before building.
5. Independent checks matter: both AeST results were re-derived by a separate reviewer before being recorded.
6. Keep failures: most of the progress came from understanding why something failed.
7. Check time dependence: a threshold that works today can be very different when structures actually formed (DBI phase change).
8. Galaxies are 'in the middle' on every simple physical axis, so no single switch can treat them differently from both smaller/earlier and larger/later structures.

## Contents

- **`camb_patches/`:** patches to [CAMB](https://github.com/cmbant/CAMB) (applied to the Sept 2026 master):
  - `camb_gboost.patch`: time-varying G in the Einstein equations;
  - `camb_varying_constants.patch`: fine-structure constant α and electron mass at recombination (`set_vconst`);
  - `camb_electron_switch.patch` (on top of the previous one): electron mass that relaxes to today's value around a chosen redshift (`set_vswitch`).
- **`cosmology_fits/`:** Cobaya fits and cosmology tests: the dark-energy law, early stretch energy, P-ACT (ACT DR6) runs, supernova shell, sky-direction and host tests, local void, H₀ without CMB, quasars, early galaxies, the Cold Spot, and a₀ from cosmic tension.
- **`galaxies_lensing/`:** SPARC rotation curves, the neighbour (external field) test, light-only predictions, KiDS lensing, lame links, mass-beats-tension, spinning capped black holes.
- **`grid_models/`:** tension-grid gravity, full tension (γ = 1), cracked mosaics (grain check), cell size vs strain, lattice gauge and fermion runs, tetrahedral grids.
- **`fluid/`:** the fluid between the layers: 1D-channel equation of state, superfluid a₀, tests of standalone superfluid dark matter.
- **`lorentz/`:** light-speed (grid spacing) limits and the AeST preferred-frame calculation.
- **`aest_upgrade/`:** our CLASS build (AeST and Khronon dark sectors + grid dark energy), fits, small-K_B instability, Khronon fits.
- **`clusters/`:** cluster tests on the Khronon base.
- **`black_holes/`:** Planck stars, lifetimes, end states, supermassive black holes.
- **`literature/`:** reviews of related papers and a survey of relativistic MOND theories.
- **`particles_collapse/`:** two electrons on a grid, the two-layer electron, acceleration- and strain-triggered collapse, flavour and neutrino studies.

## Data (not included; public sources)

- SPARC rotation curves (Lelli, McGaugh & Schombert)
- KiDS-1000 weak-lensing RAR (Brouwer et al. 2021, kids.strw.leidenuniv.nl)
- ACT DR6 CMB-only (LAMBDA)
- Planck (via Cobaya)
- DESI DR2 BAO, DES-Dovekie, Pantheon+SH0ES, Union3
- Chae et al. 2021 (ApJ 921, 104) Table 3
- Lusso et al. 2020 quasars

## Credit

- The slack rule is Milgrom's MOND (1983).
- Stored early energy belongs to the early-dark-energy family.
- The two-layer electron uses Kaplan's domain-wall construction.
- Strain-triggered collapse is Diósi–Penrose.
- A heavier early electron is a known Hubble-fix candidate.
- a₀ = cH₀/6 from the elastic response of the dark-energy medium is Verlinde's emergent gravity (2016).
- AeST is Skordis & Złośnik (2021); Khronon theory is Blanchet & Skordis (2024); the superfluid picture is Berezhiani & Khoury (2015).
- Planck stars are Rovelli & Vidotto (2014) and Haggard & Rovelli (2015).

The grid-tension framing, the dark-energy law, and the tests here are this project's.
