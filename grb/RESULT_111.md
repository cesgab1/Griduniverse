# RESULT 111 -- gamma-ray-burst distances (A118)
Files: iter111_grb.py/.txt (pre-registered form), iter111_paperform.py/.txt (the paper's own form, added after the check below).

Data check. The pre-registered expectation "b ~ 1.1-1.3, s_int 0.35-0.45" was MY slip: those numbers belong to the paper's
form (log Eiso on log Ep); my pre-registered method used the reverse (log Ep on log Eiso), where b ~ 0.5 is normal.
Run in the paper's form: gamma 1.10, s_int 0.40 dex -- as published. Data verified (plus the exact column checksums).
This adds one undeclared free choice: which way round to write the Amati relation (now 4 free choices).

| comparison (delta -2lnL vs LCDM)      | pre-registered form | paper's form |
|---------------------------------------|---------------------|--------------|
| LAW (our law)                         | -0.21               | +0.003       |
| EdS (all matter, no dark energy)      | +5.5                | +0.01        |
| BARE-OPEN (atoms only, no dark energy)| -1.1                | +9.4         |
| best Om (LCDM)                        | 0.07 (0.06-0.20)    | 0.89 (0.35-0.99) |
| Om = 0.31 penalty                     | +1.9                | +1.3         |

Verdicts against the pre-registration:
- LAW vs LCDM: indistinguishable (|delta| < 1 both ways) -- as predicted by the power check. HIT.
- Om loose and consistent with 0.31 within ~1.4 sigma both ways. HIT.
- No-dark-energy universes: the answer FLIPS with the arbitrary direction choice (EdS rejected one way, fine the other;
  BARE-OPEN favoured one way, rejected at ~3 sigma the other). So GRBs give no reliable verdict on them. Neither direction is a
  rescue; supernovae + BAO + CMB already exclude both (iterations 92-93).
Lesson: GRB Hubble diagrams depend on how the calibration relation is written -- a hidden free choice that moves Om from
0.07 to 0.89. Any GRB cosmology claim (ours included) must be shown in both forms.
