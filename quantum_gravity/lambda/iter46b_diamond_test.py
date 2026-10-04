"""
ITERATION 46b (declared fix of 46: 'nearest point in the future cone' is dominated by far points near the light cone, so it depends
on the edges of the sample -- the known infinite-valency issue -- and is not a valid frame test).
Valid test: count grid points inside causal diamonds of FIXED proper size (proper time tau between tips), oriented at rapidity
0, 1, 2 relative to the grid's frame. A grid with no preferred frame gives the same count statistics (Poisson: mean = variance
= density x volume) at every rapidity. Same expectations as 46: (S) unchanged; (T) changes with rapidity.
"""
import numpy as np
rng = np.random.default_rng(461)
rho = 400.0; T = X = 40.0
def grid(kind):
    n = rng.poisson(rho*T*X); t = rng.uniform(-T/2, T/2, n); x = rng.uniform(-X/2, X/2, n)
    if kind == "T": t = np.round(t*np.sqrt(rho))/np.sqrt(rho)
    return t, x
tau = 0.5; V = tau**2/2                                   # area of a 1+1-D diamond with proper time tau
out = ["ITERATION 46b: counts in causal diamonds of fixed proper size at different rapidities", "",
       f"expected for a frame-free grid: mean = variance = {rho*V:.1f}", ""]
for kind, name in (("S", "random in space AND time"), ("T", "random in space, synchronised ticks")):
    t, x = grid(kind)
    for eta in (0.0, 1.0, 2.0):
        c, s = np.cosh(eta), np.sinh(eta); cnt = []
        for _ in range(4000):
            t0, x0 = rng.uniform(-8, 8, 2)
            # diamond tips at +/- tau/2 along the boosted time direction (cosh, sinh)
            tl, xl = t - t0, x - x0
            tp, xp = c*tl - s*xl, c*xl - s*tl                 # coordinates in the diamond's own frame
            inside = (np.abs(xp) < tau/2 - np.abs(tp))
            cnt.append(inside.sum())
        cnt = np.array(cnt)
        out.append(f"{name:38s} rapidity {eta}: mean {cnt.mean():6.1f}   variance {cnt.var():6.1f}   variance/mean {cnt.var()/cnt.mean():.2f}")
txt = "\n".join(out); print(txt); open("iter46b_diamond_test.txt", "w").write(txt + "\n")
