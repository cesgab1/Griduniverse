"""
ITERATION 36 (QG step 3): build the conserved 'stamp' into the grid's action (unimodular-type).
Minisuperspace action of Claim 1: P = -(2/3) C (H a)^(-1/2). Promote C to a variable held constant by a Lagrange multiplier tau:
   S ⊃ -∫ dt [ C(t) * N a^3 w(a,H) - C(t) d(tau)/dt ],   w = (a H t_P)^(-1/2)
Varying tau: dC/dt = 0 (C is CONSERVED: a stamp no jostling can change). Varying C: d(tau)/dt = N a^3 w, so tau = the
w-weighted 4-volume W: C and W are a conjugate pair (as Lambda and the 4-volume are in unimodular gravity).
Consequences (expectations written before computing):
 S1 classically C is an integration constant: ANY value is allowed -> the action alone does not fix the size.
 S2 quantum uncertainty dC dW >= hbar with W = whole history gives rho ~ 1/N4 (Planck units): ~1e-245, far too small.
 S3 discreteness of cells (Poisson: dW ~ w sqrt(N4)) gives rho ~ 1/sqrt(N4): the weight w cancels -> exactly Sorkin's estimate,
    which matches the total density at whatever epoch the cells are counted (automatic, as found earlier). Counted over the
    past -> tracks the horizon (iteration 21 failure); counted over the WHOLE history -> fixed fee of the right size only if the
    whole history is ~4.4x the visible past (finite book, iteration 32).
"""
import numpy as np
N4_past = 1.9e245          # visible region, whole past (computed earlier)
obs = 1.1e-123
out = ["ITERATION 36: conserved stamp in the action (unimodular-type) -- consequences", "",
       "S1: dC/dt = 0 exactly; C is an integration constant (any value) -> no size from the classical action.",
       f"S2: hbar-uncertainty with W = 4-volume: rho ~ 1/N4 = {1/N4_past:.1e} (measured {obs:.1e}): ~1e122 x too small.",
       f"S3: cell discreteness (Poisson): rho ~ 1/sqrt(N4) = {1/np.sqrt(N4_past):.1e} for the visible past (2x measured) -- the same",
       "    automatic match; weight w cancels, so the stamp route reduces exactly to the cell tally / book.",
       "", "Reading: building the stamp into the action WORKS mathematically (C becomes conserved, so jostling cannot erase it: Coalesce's",
       "shaking point is built in) and it is ghost-free (no new propagating field). But it moves the size into the choice of",
       "which 4-volume fixes the stamp: past (fails history), whole finite history (Branch A book, untestable end; testable half =",
       "small flat wrap-around universe). The three QG routes converge on the same fork."]
txt = "\n".join(out); print(txt); open("iter36_conserved_stamp.txt", "w").write(txt + "\n")
