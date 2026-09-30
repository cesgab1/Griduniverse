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
| Early stored stretch energy | Planck + ACT DR6 | rejected by ACT (Δχ² +13 at δ = 0.03) |
| Five late-universe Hubble fixes | Pantheon+SH0ES, DESI, BBN | all ruled out (no-CMB H₀ = 68.6 ± 0.6) |

## Contents

- **`camb_patches/`:** patches to [CAMB](https://github.com/cmbant/CAMB) (applied to the Sept 2026 master):
  - `camb_gboost.patch`: time-varying G in the Einstein equations;
  - `camb_varying_constants.patch`: fine-structure constant α and electron mass at recombination (`set_vconst`);
  - `camb_electron_switch.patch` (on top of the previous one): electron mass that relaxes to today's value around a chosen redshift (`set_vswitch`).
- **`cosmology_fits/`:** Cobaya fits and cosmology tests: the dark-energy law, early stretch energy, P-ACT (ACT DR6) runs, supernova shell, sky-direction and host tests, local void, H₀ without CMB, quasars, early galaxies, the Cold Spot, and a₀ from cosmic tension.
- **`galaxies_lensing/`:** SPARC rotation curves, the neighbour (external field) test, light-only predictions, KiDS lensing, lame links, mass-beats-tension, spinning capped black holes.
- **`grid_models/`:** tension-grid gravity, full tension (γ = 1), cracked mosaics (grain check), cell size vs strain, lattice gauge and fermion runs, tetrahedral grids.
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

The grid-tension framing, the dark-energy law, and the tests here are this project's.
