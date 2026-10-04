"""
Does the wobble in the leftover dark energy line up with a change in HOW it was measured? DESI uses different kinds of galaxies
at different redshifts: z 0.51, 0.71 = luminous red galaxies (LRG); 0.93 = LRG + emission-line galaxies combined; 1.32 = ELG;
1.48 = quasars; 2.33 = Lyman-alpha forest. Compare the leftover from LRG-only bins with the rest (full covariance by Monte Carlo).
"""
import numpy as np
src = open("leftover.py").read().split("# Claude 1 prediction")[0].split("# Claim 1 prediction")[0]
g = {}; exec(src, g); X, z = g["X"], g["z"]
lrg = np.isin(np.round(z, 2), [0.51, 0.71]); rest = ~lrg
def wmean(cols):
    Xs = X[:, cols]; C = np.cov(Xs.T); w = np.linalg.inv(np.atleast_2d(C)); one = np.ones(cols.sum())
    return (one @ w @ Xs.T)/(one @ w @ one)          # per-sample weighted mean
A, B = wmean(lrg), wmean(rest); D = A - B
out = [f"LRG bins (z 0.51, 0.71): leftover {A.mean():.3f} +/- {A.std():.3f}",
       f"other tracers (z 0.93-2.33): leftover {B.mean():.3f} +/- {B.std():.3f}",
       f"difference {D.mean():+.3f} +/- {D.std():.3f}  ({D.mean()/D.std():.1f} sigma)",
       "The step sits exactly where the galaxy type changes (LRG -> LRG+ELG at z ~ 0.9), and also near where Claim 1 puts the",
       "dark-energy peak (z ~ 0.68) and the start of acceleration. Data alone cannot separate 'tracer effect' from 'physics'."]
txt = "\n".join(out); print(txt); open("tracer_split.txt", "w").write(txt + "\n")
