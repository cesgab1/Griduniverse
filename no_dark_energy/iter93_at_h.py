# BAO and CMB-angle pulls of the bare universes at the measured expansion rates (h not fitted)
import os, sys; sys.argv = ["x"]
src = open("iter93_bare.py").read().split('out = [f"Iteration 93')[0]; exec(src)
for m in ["BARE-OPEN", "BARE-FLAT"]:
    for h in (0.674, 0.73):
        o = build(m, h); cb, cs, lp = chis(o)
        print(f"{m:9s} h={h}: BAO chi2 {cb:9.1f} (13 points), r_d {o['rd']*RD_CAL:.1f} Mpc, l_A {o['lA']:.1f} ({lp:+.0f} sigma)")
