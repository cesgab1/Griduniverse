# PREREG 104 -- do the largest-scale anomalies point toward one direction (a far-off centre)? (Coalesce)
Committed BEFORE running iter104_align.py. NOT BLIND: I roughly know these directions already.
Source: Aluri et al. 2023 (CQG 40 094001, arXiv:2207.05765) Table 1, spot-checked (CMB dipole, CatWISE, power asymmetry match).
One entry per independent phenomenon (duplicates and known mass concentrations removed):
| key | phenomenon | (l, b) | type |
| CMB | CMB dipole (usually our motion) | (264.0, 48.3) | vector |
| RADIO | radio-galaxy count dipole NVSS | (248, 44) | vector |
| QSO | quasar count dipole CatWISE | (238.2, 28.8) | vector |
| LOWL | CMB octopole axis (quadrupole aligned) | (239, 64.3) | axis |
| HPA | CMB hemispherical power asymmetry (Planck DM) | (227, -15) | vector |
| PARITY | CMB odd mirror-parity maximum | (264, -17) | axis |
| CLUSTER | galaxy-cluster H0 anisotropy (Migkas) | (280, -15) | vector |
| BULK | bulk flow (Atrio-Barandela) | (282, 6) | vector |
| SN | our supernova near-shell dipole (iteration 103) | (316.5, 12.6) | vector |
| PARAM | cosmological-parameter dipole (Yeung) | (48.8, -5.6) | vector |
| COLD | CMB cold spot (from memory) | (209, -57) | vector |
Statistic: length R of the average unit vector (axes take the sign that best aligns). p-value from 100,000 random isotropic
sets with identical rules. Version 1: vectors as published. Version 2: everything sign-free (axes).
Also: best common direction; which entries are within 45 deg of it.
Expectations: version 2 p < 0.01 (most sit in l 210-320); version 1 p 0.01-0.1. A far-off centre predicts tight alignment of ALL;
shared alignment near the CMB dipole could equally come from our own motion or survey systematics -> cannot prove a centre.
Look-elsewhere: list compiled by an anisotropy review; selection effects inflate significance.
