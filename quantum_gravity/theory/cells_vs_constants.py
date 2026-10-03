"""
If the Planck-size cells themselves stretched (cell size ∝ a^s_cell), the Planck length l_P^2 = hbar G / c^3 would grow.
How fast may it grow? Use only G (hbar, c held fixed; any split makes at least one constant vary as much).
"""
import numpy as np
H0_yr = 67.5/3.0857e19*3.156e7         # H0 in 1/yr
# Lunar laser ranging: Gdot/G = (7.1 ± 7.6)e-14 /yr  -> 2-sigma |Gdot/G| < 2.23e-13 /yr
lim_llr = (7.1 + 2*7.6)*1e-14
s_llr = lim_llr/(2*H0_yr)
# Big-bang nucleosynthesis: |dG/G| < ~10% at z ~ 4e8 ;  l_P^2 ∝ a^(2 s) -> G(BBN)/G0 = (1+z)^(-2 s)
s_bbn = np.log(1.1)/(2*np.log(4e8))
# our CMB constant test: alpha_rec/alpha0 = 1.0031 ± 0.0013 (alpha ∝ 1/c); if c ∝ cell size ∝ a^s -> c_rec/c0 = 1101^-s
s_cmb = np.log(1.0031 + 2*0.0013)/np.log(1101)
out = [f"cell-stretch exponent allowed by lunar laser ranging (G today):      s_cell < {s_llr:.1e}",
       f"cell-stretch exponent allowed by nucleosynthesis (G in the first minutes): s_cell < {s_bbn:.1e}",
       f"cell-stretch exponent allowed by the CMB (c, via alpha, if c = cell/tick): s_cell < {s_cmb:.1e}",
       "dark-energy law fit (stretch_scan): s_DE = 0.91-1.03 ± 0.15-0.20",
       "=> the length that stretches in the dark-energy law is NOT the Planck cell (which must stay fixed to <1e-3 in exponent,",
       "   i.e. cells ARE being added). The counting must run over a COMOVING length, distinct from the cell size."]
txt = "\n".join(out); print(txt); open("cells_vs_constants.txt", "w").write(txt + "\n")
