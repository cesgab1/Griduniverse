# The size of dark energy

**Iteration 21: random spin flips of grid cells, with the grid's memory** (iter21_spin_nudges.py/.txt)

**Size (good).**
- The flips remembered now number N = 6×10²⁴¹. Their leftover imbalance gives 75× the critical density for an order-one
  coefficient.
- Observed: 0.69×. So the size problem's 120 orders of magnitude shrink to a factor of about 100.
- This is Sorkin's everpresent-Λ scaling with the toy's memory window (BORROWED).

**History (fails).**
- A size recomputed from the remembered region at every epoch TRACKS the critical density, so dark energy becomes a
  constant fraction of the universe (Ω_DE = c × random).
- The CMB needs early dark energy ≲ 0.02 near recombination, and nucleosynthesis ≲ 0.1.
- Monte Carlo of 20,000 histories: no amplitude passes the early limits AND gives Ω_DE = 0.62–0.76 today. The pass rate
  for both is 0.0000%.
- Even passing amplitudes would give a noisy w(z) (a jump every ~⅓ e-fold), unlike the smooth single crossing the data
  prefer. Half of all draws are negative.

**Lesson.** Size and shape are linked. A size set by "what fits in our horizon now" makes dark energy track everything
else. Our law needs it nearly constant while matter thins out. So the size must be ACCUMULATED or FROZEN IN, not
recomputed from the horizon each moment.

**Candidates (OPEN):**
- **Frozen at launch:** the leftover imbalance is set once, when the quiet grid started stretching (iteration 20):
  ρ_DE ~ ρ_P/√N_launch, constant afterwards, with its shape modulated by the jostling. It needs N_launch ~ 10²⁴⁴ cells:
  a large quiet grid before the Big Bang (emergent-universe picture). It loses the automatic "why now".
- **Only flips that are never undone count** (a ratchet). This needs a rule for which flips persist.

Each must be computed, not tuned, under the self-audit rules.

**Iteration 22: built-in tension with a yield threshold (stick-slip; Coalesce's rubber)** (iter22_stick_slip.py/.txt)

**Setup.** Each cell's strain grows at the stretch rate. At the threshold ε_c a new cell is added and the strain drops by
ε_c. Simulated over a matter → Λ history with 20,000 cells.

**Results.**
- **Prompt relief:** the leftover strain stays CONSTANT (⟨ε²⟩/(ε_c²/3) = 0.985 / 0.999 / 0.998 / 1.001 at
  z = 29 / 3 / 1 / 0). A frozen size comes naturally, and it behaves exactly like a cosmological constant (w = −1).
- **Short relief delay:** an overshoot that is bigger early (1.5–2.9× at z = 29) and fades as the expansion slows. So
  w > −1 always and never crosses −1, against the data's preferred crossing.
- **Long delay:** relief cannot keep up and strain piles up without limit (runaway). Not a viable grid.
- **Size:** needs ε_c = 8×10⁻⁶², i.e. one cell of stretch per row of 1.2×10⁶¹ cells (1.9×10²⁶ m), which is 1.4× today's
  Hubble length.
  - If that row follows the current horizon, the size tracks H² (iteration 21's failure).
  - If it was fixed once, its length equals today's horizon by coincidence ("why now").

**Data.**
- Prompt stick-slip = Λ, so it loses Claim 1's Δχ² gain (−5.4 / −7.2 / −6.8).
- As a positive baseline under the jostled toy it is disfavoured (iteration 13).

**Lesson.**
- Built-in tension does give a FROZEN size, which fixes iteration 21's history problem. Coalesce's point holds.
- But its natural behaviour is a pure constant. It does not carry Claim 1's shape, and the small threshold is the old
  size problem in a new form.
- Across 21 and 22: every natural size mechanism ties the size to the horizon (right size now, wrong history) or freezes
  it at a value that must match today's horizon by coincidence.

**Iteration 23: size tied to the horizon (Coalesce), via holographic dark energy (Li 2004, future event horizon)**
(iter23_holographic.py, iter23_raw.txt)

Fit to the same data as Claim 1, with one extra parameter c. Δχ² vs Λ:

| supernovae | holographic (1 extra parameter) | Claim 1 (0 extra) |
|---|---|---|
| Pantheon+ | +44.5 | −5.4 |
| DES-Dovekie | +45.5 | −7.2 |
| Union3 | +41.5 | −6.8 |

- Best c = 0.69–0.73, which gives w today = −1.10 to −1.14: phantom NOW, and crossing the opposite way to the data.
- All three datasets object (DES-Dovekie case vs Λ: BAO +16, supernovae +5, CMB +24).
- A first run that switched dark energy off above z = 30 was buggy (iter23_raw_TRUNCATED_BUG.txt), because holographic dark
  energy is not negligible early. The full-history rerun changes nothing material.

**Lesson.** Tying the size to the future horizon, the most established version of Coalesce's idea, gets the order of
magnitude right but the history badly wrong (about 50 in χ² worse than Claim 1). The data strongly favour Claim 1's order:
phantom in the past, w > −1 now. A horizon-based size would need a different shape mechanism (such as the grid's memory),
not the holographic one.

## Iteration 24: the trampoline (Pools sag the sheets), PRE-REGISTERED before fitting
Motivation from the sky: dark matter is constant to a few % since the CMB, so dark energy cannot be drained from the Ocean;
only the Ocean's ARRANGEMENT (pooling into halos) can change. Predicted histories (`iter24_trampoline.py`,
`iter24_part1_predictions.txt`, `iter24_predictions.json`), written before any data comparison:
- T1 (tension ~ fraction of Ocean in Pools, M_min 1e10): dark energy only GROWS -> w < -1 at all z, NO crossing,
  w0 = -1.09, bins -1.14 ... -1.44. Wrong direction vs the DR2-era hint (w > -1 today).
- T2 (tension ~ total sag energy rho_m <|Phi|>): peaks at z ~ 2, fades since -> w0 = -0.19, crossing z ~ 2.0. Far too fast.
- Expectation written now: both worse than Lambda; T2 badly. Bare size: sag = 4.6e-7 of the need (amplifier ~2e6, not derived).
- The galaxy-scale (local) version of the trampoline is Family 1 (static links -> plain Newton); not re-run.

### Iteration 24 results (part 2, after the commit 2abb7e2 above)
Fit to DESI DR2 BAO + Planck distance priors + each supernova set; fixed shape, no dark-energy parameter
(`iter24_part2_fits.txt`; new `TAB` option in fit_law.py, validated: a flat shape reproduces Lambda to 0.00).

| | Pantheon+ | DES-Dovekie | Union3 |
|---|---|---|---|
| T1 (fraction in Pools) | +41.3 | +50.4 | +31.5 |
| T2 (sag energy) | +513 | +596 | +282 |

(Delta chi^2 vs Lambda; Claim 1 is -5.4 / -7.2 / -6.8.) Both excluded, as pre-registered. Bare size: 4.6e-7 of the need.

**Lesson (the useful part).** Any tension set by clumping follows the STRUCTURE clock. Clumping only grows (T1: dark energy
always rising, w < -1, the opposite of the hint), or its energy peaks at z ~ 2 and drains as matter thins (T2). The data want
dark energy that peaked near z ~ 0.5-0.7 and is fading: Claim 1's peak sits exactly at the onset of acceleration (q = 0).
So the clock the data read is the EXPANSION clock (adot), not the structure clock. Dark matter enters only through its
thinning (it sets when q crosses 0), not through its pooling. A size mechanism must be tied to the expansion rate, yet not
track H^2 (iteration 21). Mixing T1 and T2 to place a peak at z ~ 0.6 would be an accommodation (new free weight, no new
prediction), so it was not done.
