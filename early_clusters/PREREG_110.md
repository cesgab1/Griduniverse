# PREREG 110 -- population count: early ultra-massive galaxies vs standard cosmology (Coalesce: count populations)
Committed BEFORE running iter110_population.py.
Data (Xiao et al. 2024, Nature 635, 311; arXiv:2309.02492, verified): FRESCO, 124 arcmin^2 (GOODS-N + GOODS-S), z = 4.9-6.6,
spectroscopic. Three galaxies: log M* = 11.37 (z 5.58), 11.18 (z 5.31), 11.04 (z 5.18). Observed count with log M* >= 11.04: 3.
Method: our Planck-cosmology halo abundance (Sheth-Tormen, iteration 109 code). A galaxy of stellar mass M* needs a halo of mass
>= M* / (eps x f_b), f_b = 0.157, eps = share of the halo's ordinary matter turned into stars. Expected number in the FRESCO volume
for eps = 1 (absolute maximum), 0.5, 0.2 (typical upper value at later times); Poisson chance of >= 3.
Sensitivity: masses 0.15 dex lower (typical error scatter pushing estimates up, 'Eddington bias').
Plus JADES-ID1 (iteration 109) as a separate object class in an overlapping field.
Expectations:
- eps = 1: expected >~ 3 -> standard cosmology itself is NOT broken.
- eps = 0.2: expected ~0.01-0.1 -> chance of >= 3 below 1e-4 -> galaxy formation must be ~2-5x more efficient early on.
- Populations do not rescue a 'fluke' reading: the tension is with how galaxies form, not with gravity or dark matter's amount.
