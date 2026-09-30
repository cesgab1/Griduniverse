import sys; sys.argv = ["x", sys.argv[1] if len(sys.argv) > 1 else "0.315"]
src = open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/void_test.py").read()
src = src.split("    res = {}")[0]   # keep setup + functions for the first redshift column only
src = src.replace('for zcol in ("zCMB", "zHD"):', 'for zcol in ("zCMB",):')
exec(src)
from scipy.optimize import minimize_scalar
bp = minimize_scalar(lambda H: chi2(H, 0, 1), bounds=(60, 80), method="bounded")
print(f"Om={Om}  B' no void, H0 with ruler prior: chi2 = {bp.fun:.1f} (H0 {bp.x:.2f})   [A, no prior: {minimize_scalar(lambda H: chi2(H,0,1,False), bounds=(60,80), method='bounded').fun:.1f}]")
for R in (300, 500, 800):
    for eps in (0.04, 0.075):
        f = minimize_scalar(lambda H: chi2(H, eps, R), bounds=(60, 80), method="bounded")
        print(f"   void R={R} Mpc forced local excess {eps*100:.1f}%: chi2 = {f.fun:.1f} (H_out {f.x:.2f})  vs no-void-with-prior {bp.fun:.1f}: {f.fun-bp.fun:+.1f}")
for R in (800, 1600, 2400, 4000):
    b = minimize(lambda p: chi2(p[0], p[1], R), [67.9, 0.03], method="Nelder-Mead")
    print(f"   void R={R}: best excess {b.x[1]*100:+.1f}%, H_out {b.x[0]:.2f}, chi2 {b.fun:.1f} ({b.fun-bp.fun:+.1f})")
