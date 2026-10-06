# RESULT 106 -- is the hemisphere split real? (survey drop-one + independent Union3 data)
Output: iter106_checks.txt
Check 1 -- Pantheon+, drop one large survey at a time (z > 0.05; coarser fit grid and a new set of 300 random axes):
| sample | toward/away SNe | matter-share difference (toward - away) | random-axis p |
| all surveys | 271/785 | -0.060 +/- 0.032 | 8.3% |
| drop survey 1 | 252/483 | -0.070 +/- 0.039 | 8.7% |
| drop survey 4 | 233/663 | -0.070 +/- 0.042 | 11.3% |
| drop survey 10 | 207/646 | -0.070 +/- 0.039 | 8.7% |
| drop survey 15 | 158/632 | -0.050 +/- 0.039 | 12.3% |
The split does not hinge on any single survey (EXPECTATION MISSED: I predicted one survey would drive it). But the same baseline
gives p = 8.3% here vs 3.7% in iteration 105 (finer grid, 1500 axes): the significance is fragile to analysis details.
Check 2 -- Union3 inputs, independent light-curve calibration (simplified standardisation, 1497 SNe at z > 0.05):
dark-energy share toward 0.530, away 0.570 -> OPPOSITE sense to Pantheon+, difference +0.040 +/- 0.056 in matter share, p = 64%.
(Absolute values are off from usual because the method is simplified; the toward/away comparison uses the same method on both sides.)
Expectation: 'same sense, weaker, p > 0.1' -- sense MISSED (opposite), p > 0.1 MET.
VERDICT: the hint of stronger speed-up toward (260, +12) is NOT reproduced by an independently calibrated catalogue and its strength
in Pantheon+ depends on analysis details. Best reading: a statistical fluctuation (or Pantheon+-specific), not a direction-dependent
dark energy. The far-off-centre picture gains no support from supernovae; it is not ruled out (centre could be farther).
