import numpy as np
from scipy.stats import poisson, norm
exec(open("iter109_rarity.py").read().split("out = []")[0])
fb = 0.157
obs = [11.37, 11.18, 11.04]; Mmin_star = 10**min(obs)
out = [f"FRESCO 124 arcmin^2, 4.9 < z < 6.6; observed galaxies with log M* >= {np.log10(Mmin_star):.2f}: {len(obs)}"]
for shift, lab in ((0.0, "masses as published"), (-0.15, "masses 0.15 dex lower (error-scatter correction)")):
    out.append(f"-- {lab}")
    for eps in (1.0, 0.5, 0.2):
        Mh = Mmin_star*10**shift/(eps*fb)
        lam = N_above(Mh*h, 4.9, 6.6, 124)
        p = poisson.sf(len(obs) - 1, lam)
        sig = norm.isf(p) if p > 0 else np.inf
        out.append(f"   eps = {eps:.1f}: needs halo >= {Mh:.1e} Msun; expected {lam:.2e}; chance of >= 3: {p:.1e} ({sig:.1f} sigma)")
txt = "\n".join(out); print(txt); open("iter110_population.txt", "w").write(txt + "\n")
