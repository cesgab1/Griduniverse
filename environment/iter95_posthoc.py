"""POST-HOC (not pre-registered) robustness check of part B: medians and 3-sigma-clipped means."""
import numpy as np
exec(open("iter95_environment.py").read().split("# C. Milky Way")[0].replace("print(txt)", "pass"))
def clip(x):
    m = np.ones(len(x), bool)
    for _ in range(5):
        mu, s = x[m].mean(), x[m].std(ddof=1); m = np.abs(x - mu) < 3*s
    return x[m]
uc, fc = clip(u), clip(f)
d2 = uc.mean() - fc.mean(); e2 = np.hypot(uc.std(ddof=1)/np.sqrt(len(uc)), fc.std(ddof=1)/np.sqrt(len(fc)))
out = [f"POST-HOC: medians  UMa {np.median(u):+.3f}  field {np.median(f):+.3f}  diff {np.median(u)-np.median(f):+.3f} dex",
       f"POST-HOC: 3-sigma clipped  UMa N={len(uc)} {uc.mean():+.3f}  field N={len(fc)} {fc.mean():+.3f} (scatter {fc.std(ddof=1):.3f})  diff {d2:+.3f} +/- {e2:.3f} ({d2/e2:+.1f} sigma)"]
print("\n".join(out)); open("iter95_posthoc.txt", "w").write("\n".join(out) + "\n")
