# Same pre-registered comparison, measured: a0 of the deepest 10% of SPARC galaxies vs those at median depth (40-60%).
import numpy as np
exec(open("../missing_input/iter96_missing_input.py").read().split("names = dict")[0])
x = np.log10([r_["v2"] for r_ in rows])
mid = (x > np.percentile(x, 40)) & (x < np.percentile(x, 60)); deep = x >= np.percentile(x, 90)
rng = np.random.default_rng(2)
def med(a): return np.median(a)
r = 10**(med(la[deep]) - med(la[mid]))
bs = [10**(med(rng.choice(la[deep], deep.sum())) - med(rng.choice(la[mid], mid.sum()))) for _ in range(2000)]
lo, hi = np.percentile(bs, [16, 84])
pred = 10**(0.327*(np.median(x[deep]) - np.median(x[mid])))
out = f"deepest 10% (N={deep.sum()}) vs median-depth (N={mid.sum()}): measured boost {r:.2f} [68%: {lo:.2f}-{hi:.2f}]; depth law predicts {pred:.2f}"
print(out); open("iter101_groups.txt", "a").write(out + "\n")
