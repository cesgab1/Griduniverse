"""
ITERATION 28: can the number of sheets (layers) tell how long the 'book' (total pages x) runs? (Coalesce, Oct 2026)
Each cell = coordinates + time stamp (layer, page). Layer window from earlier work: 1e7 <~ N <~ 1.2e15.
Two physically motivated links, stated BEFORE computing (no other rules tried, to avoid number-hunting):
 M1 sheets ARE pages: x = N, each page one tick of the grid (Planck time, or the cell time sqrt(N) t_P).
 M2 sheets set the tick: with N layers the effective cell is sqrt(N) l_P (species bound), its energy density rho_P/N^2,
    and the book size rule rho_DE = rho_cell / x^2 with x counted in cell ticks.
PASS criterion: the implied total lifetime must exceed today's age (13.8 Gyr) and give the measured dark-energy size.
"""
import numpy as np
tP = 5.391e-44; t0 = 13.8e9*3.156e7; rhoDE_over_rhoP = 1.1e-123      # dark energy in Planck units (law_from_grid.md)
n_now = t0/tP
x_needed = np.sqrt(1/rhoDE_over_rhoP)
out = ["ITERATION 28: sheets -> book length", "",
       f"pages read so far (Planck ticks): {n_now:.2e}; book length needed for the dark-energy size (rule 1/x^2): {x_needed:.2e}",
       f"  -> total lifetime {x_needed*tP/3.156e16:.0f} Gyr, i.e. {x_needed/n_now:.1f} x today's age", ""]
for N in (1e7, 3e10, 1.2e15):
    T1a = N*tP; T1b = N*np.sqrt(N)*tP
    T2 = np.sqrt(1/(N*rhoDE_over_rhoP))*tP          # rho_P/N^2 * (sqrt(N) tP / T)^2 = rho_DE
    out.append(f"N = {N:.1e}:  M1 book lasts {T1a:.1e} s (Planck ticks) or {T1b:.1e} s (cell ticks);  "
               f"M2 book lasts {T2/3.156e16:.2e} Gyr = {T2/t0:.1e} x today's age")
out += ["", "Verdict: M1 fails by > 40 orders of magnitude (the book would have ended a fraction of a second after the start).",
        "M2 fails too: more sheets make the book SHORTER (lifetime ~ 1/sqrt(N)); even N = 1e7 ends ~800x too soon.",
        "Reading: the sheets exist side by side, all at once (they average the jostling NOW); they are not pages in sequence.",
        "The book's length is a separate number from the sheet count. A fixed total CELL count (space x layers x pages) is the",
        "grid version of sequestering's 'total spacetime volume', but it also needs the size of space beyond our horizon (unknown)."]
txt = "\n".join(out); print(txt); open("iter28_sheets_to_pages.txt", "w").write(txt + "\n")
