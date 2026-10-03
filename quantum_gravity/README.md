# Quantum gravity and the Grid Universe: gap map + one test

`area_law.py` → `area_law.txt`: each node is a quantum oscillator, coupled through the tension links (w_ij = A_ij/d_ij).
In the ground state, the entanglement entropy of a ball of cells scales as R^2.00. It tracks the wall area cut by the
boundary, as black-hole entropy does, not the volume.
The coefficient (0.060 per cell-wall area) is not the Bekenstein-Hawking 1/4. It depends on cell size and on how many
fields live on the grid, so it is not a prediction. This reproduces Bombelli et al. 1986 / Srednicki 1993 on our random
mosaic (BORROWED).

The gap map is in IDEAS_LEDGER.md ("Quantum gravity: where the model's gaps are").

## Graviton tests (gap 2: can the grid carry a spin-2 wave?)

**Counting.** One number per link (tension) carries helicity 0 only. A displacement per node (springs) carries 0 and ±1.
Only a tensor per cell carries ±2 (the graviton). A tensor per cell is the cell's SHAPE, which is fixed by its link
lengths. The graviton is a ripple in cell shapes: the + polarisation stretches cells one way and squeezes them the other.

**`graviton_test.py` (derivative method): FAILS on the random grid.** Gauge residual 0.9 and growing toward long
wavelengths; the cubic lattice is exact (1e-16). The failure is in the method: derivatives of derivatives amplify the
grid's cell-to-cell irregularity.

**`regge_gauge.py` (geometric method: link lengths as the variables, Regge calculus 1961): EXACT on the random grid.**
The action is Σ (link length × deficit angle), and curvature is the angle missing around a link. The Hessian is computed
by complex-step differentiation.

| grid | zero-cost directions found / expected (3 × interior nodes) | gauge residual vs random |
|---|---|---|
| fully random | 375 / 375 | 2.5e-13 vs 3.2e-2 |
| shaken lattice | 511 / 510 | 2.9e-13 vs 2.7e-2 |

So "moving the nodes" (relabelling) costs nothing, exactly, however irregular the grid. This is the symmetry that leaves
gravity only its 2 polarisations. In 3-D there are no gravitational waves; the 4-D graviton of linearised Regge calculus
is BORROWED (Rocek & Williams 1981).

**Lesson for the picture.** Moving nodes is not gravity; it is relabelling. Gravity is when the link lengths can no
longer be fitted by any placement of the nodes in flat space (deficit angle ≠ 0). The tension-link equation (E1) is only
the Newtonian shadow of this.

**Prior art for the quantum version.** A random simplicial grid plus a preferred time slicing (our stamp), summed over
quantum mechanically, is Causal Dynamical Triangulations (Ambjørn, Jurkiewicz & Loll 2004-08: a 4-D de Sitter-like
universe emerges; dimension runs to ~2 at short distances). Its continuum cousin is Hořava gravity, whose low-energy
limit is the khronometric theory our Khronon base belongs to. Our quantum gravity is therefore BORROWED from CDT/Hořava.
Open question that could be ours: does a CDT/Hořava-type grid produce the Tension-Rate Law (ρ_DE ∝ ȧ^(−1/2))?
