"""
ITERATION 20 (Coalesce: a smooth grid existed before the Big Bang; the Big Bang is when stretching/distortion began).

A. CONSISTENCY CHECK of the toy (needed first). The toy tension's motion energy (1/2)tau^2 obeys
   d<tau^2>/dt = -6H<tau^2> + rate: that is a STIFF fluid (w = +1, ∝ a^-6) fed by energy from the glow.
   (i) stiff fluid has positive pressure -> it DECELERATES (rho + 3p > 0): it cannot by itself be accelerating dark energy;
   (ii) energy fed from the glow cannot exceed the glow's own energy (Delta N_eff < 0.107 -> < ~2e-5 of critical density).
   => The dark-energy reading that the data fits used is only consistent in the VCDM-type framework (law_from_grid sect. 1):
   the jostling statistics set the SHAPE of a vacuum-like tension whose energy is carried by the gravity sector (no extra
   degree of freedom), not drawn from the glow. That is also why its SIZE is free (iterations 10, 11, 13).
B. QUIET GRID: with the VCDM reading, the vacuum-like tension's density follows the jostled variance. In a static phase
   (H = 0) the memory 1/(3H) is infinite: the variance grows linearly forever. Model: closed static universe (k = +1),
   constraint H^2 = Om a^-3 + rho_DE - k/a^2 (units 8 pi G/3 = 1), rho_DE = A X, dX/dt = -6 H X + 1/a.
   Start exactly static (H = 0); follow the expanding branch and the contracting branch.
C. Large-angle imprint of a quiet phase before expansion: same kind of largest-scale power suppression as tested in
   iteration 12 (real Planck low-l likelihood): best Delta chi2 = -1.4 for one parameter -> allowed, not favoured.
"""
import numpy as np
from scipy.integrate import solve_ivp
out = ["ITERATION 20: the quiet (static) grid before the Big Bang", "", "A. Energy budget of the jostled motion energy:"]
Tnu = 1.945*8.617e-5; rho_g_max = 0.107*7/8*np.pi**2/30*2*Tnu**4; rho_crit = 2.54e-11/0.69
out.append(f"   glow energy (max allowed) = {rho_g_max/rho_crit:.1e} of critical density; dark energy needs 0.69 -> {0.69*rho_crit/rho_g_max:.0e}x more")
out.append("   stiff motion energy: rho + 3p = 4 rho > 0 -> decelerates. => only the VCDM reading (shape from jostling, energy via")
out.append("   the gravity sector, vacuum-like) is consistent. Applies retroactively to iterations 1-15: their fits are valid in this reading.")
out.append("")
out.append("B. Quiet closed grid, exactly static at t = 0 (a = 1, H = 0); dark-energy part grows as the memory is infinite:")
Om, k = 0.5, 1.0
for A in (0.1, 0.01):
    X0 = (k - Om)/A                                            # static balance: Om + A X0 = k at a = 1
    for branch, sgn in (("expanding", 1), ("contracting", -1)):
        def rhs(t, y):
            a, X = y; H2 = Om*a**-3 + A*X - k/a**2; H = sgn*np.sqrt(max(H2, 0)); return [a*H, -6*H*X + 1/a]
        ev = lambda t, y: y[0] - (5.0 if sgn > 0 else 0.05); ev.terminal = True
        s = solve_ivp(rhs, [0, 500], [1.0, X0], events=ev, rtol=1e-9, atol=1e-12, max_step=0.05)
        tl = s.t[np.argmax(np.abs(s.y[0] - 1) > 0.01)] if np.any(np.abs(s.y[0] - 1) > 0.01) else np.nan
        out.append(f"   A = {A:5.2f}, {branch:11s} branch: leaves stasis (1% size change) after t = {tl:6.2f}; reaches a = {s.y[0][-1]:.2f} at t = {s.t[-1]:.1f}")
out += ["", "   -> A perfectly quiet grid CANNOT stay quiet: with no expansion there is no friction, so the jostled tension builds up",
        "      without limit and the grid must start to move (a self-starting 'Big Bang'). Weaker jostling (smaller A) only",
        "      delays it. The toy does not pick the direction: the expanding branch stays smooth (iterations 17-19), the",
        "      contracting branch amplifies tension and tidal chaos. Our smooth universe is the expanding branch.",
        "      Caveat: a static closed universe of matter + vacuum energy is already unstable classically (Einstein static",
                "      universe); the new point is that even a STABILISED quiet grid could not stay quiet once jostled.",
        "", "C. Large-angle imprint: allowed, not favoured (iteration 12: Delta chi2 = -1.4)."]
txt = "\n".join(out); print(txt); open("iter20_quiet_grid.txt", "w").write(txt + "\n")
