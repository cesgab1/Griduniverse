"""POST-HOC robustness for iteration 96 (added after seeing a non-zero least-squares slope): binned medians, Theil-Sen slope, bootstrap."""
import numpy as np
from scipy.stats import theilslopes
exec(open("iter96_missing_input.py").read().split("names = dict")[0])
rng = np.random.default_rng(1)
out = []
for k, lab in (("M", "baryonic mass"), ("v2", "well depth v^2"), ("R", "size")):
    x = np.log10([r_[k] for r_ in rows]); xc = np.log10(cl[k])
    q = np.percentile(x, [0, 25, 50, 75, 100])
    meds = [np.median(la[(x >= a) & (x <= b)]) for a, b in zip(q[:-1], q[1:])]
    ts = theilslopes(la, x)[0]
    bs = [theilslopes(la[i], x[i])[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(300))]
    e = np.std(bs); dx = xc - np.median(x)
    nr = [np.log10(t)/dx for t in (5.3, 17)]
    out.append(f"{lab:15s} quartile medians a0: " + " ".join(f"{10**m:.2e}" for m in meds))
    out.append(f"{'':15s} Theil-Sen slope {ts:+.3f} +/- {e:.3f} (bootstrap); needed {nr[0]:+.2f} (x5.3) / {nr[1]:+.2f} (x17)"
               f" -> tension {(nr[0]-ts)/e:+.1f} / {(nr[1]-ts)/e:+.1f} sigma")
    out.append(f"{'':15s} extrapolated cluster a0 = {10**(np.median(la) + ts*dx):.2e} m/s^2 (x{10**(ts*dx):.1f})")
txt = "\n".join(out); print(txt); open("iter96_robust.txt", "w").write(txt + "\n")
# POST-HOC 2: is the trend a smooth power law, or only the smallest galaxies? Theil-Sen on the upper 75% only.
out2 = []
for k, lab in (("M", "baryonic mass"), ("v2", "well depth v^2"), ("R", "size")):
    x = np.log10([r_[k] for r_ in rows]); m = x > np.percentile(x, 25)
    ts = theilslopes(la[m], x[m])[0]
    bs = [theilslopes(la[m][i], x[m][i])[0] for i in (rng.integers(0, m.sum(), m.sum()) for _ in range(300))]
    dx = np.log10(cl[k]) - np.median(x[m])
    out2.append(f"{lab:15s} upper 75% only: slope {ts:+.3f} +/- {np.std(bs):.3f} -> cluster a0 x{10**(ts*dx):.1f} (needed x5.3-17)")
print("\n".join(out2)); open("iter96_robust.txt", "a").write("\n" + "\n".join(out2) + "\n")
