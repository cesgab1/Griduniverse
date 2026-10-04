# ITERATION 52b: does the eps-dependence of the curvature reader converge with resolution? (S3, eps = 0.2 vs 0.6)
import numpy as np, sys, time
sys.argv = ["x", "5"]
src = open("iter52_curvature_squared.py").read().split("m = int(sys.argv[1])")[0]
exec(src)
for m in (5, 7, 9, 11, 13):   # m = 15 ran out of memory (8 GB)
    t0 = time.time(); R = {}
    for e in (0.2, 0.6):
        R[e] = 30*kappa("S3", None, e, m)/e**2
    print(f"m={m}: R(0.2)={R[0.2]:.5f} R(0.6)={R[0.6]:.5f} ratio R(0.6)/R(0.2)-1 = {R[0.6]/R[0.2]-1:+.5f}  ({time.time()-t0:.0f}s)", flush=True)
