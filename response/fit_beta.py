import os, json, numpy as np
from scipy.optimize import minimize
os.chdir("/home/claude/griduniverse/cosmology_fits")
exec(open("fit_law.py").read().split("Ndata = len(bao)")[0])
out = {}
grid = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5]
prof = []
for b in grid:
    BETA = b
    r = best("LAW", [[0.31, 0.68, 0.0224], [0.30, 0.67, 0.0223]])
    prof.append((b, float(r.fun), [float(v) for v in r.x]))
    print("beta", b, r.fun, r.x, flush=True)
print("RESULT", SNSET, json.dumps(prof))
