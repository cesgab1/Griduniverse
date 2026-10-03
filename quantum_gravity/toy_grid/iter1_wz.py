"""Where does w cross -1 when the tension has memory 1/(kappa H)? (Om 0.31, h 0.68; lin = string-like, quad = spring-like)"""
import numpy as np, os
os.environ["SNSET"] = "UNION3"; import sys
src = open("iter1_memory.py").read().split("# ---------- Part A")[0]
g = {"__name__": "x", "__file__": os.path.abspath("iter1_memory.py")}; exec(src, g)
a = np.linspace(1/31, 1, 4000); out = ["memory test: w(z) of the grid tension, Om = 0.31, h = 0.68", " energy kappa   w0     z_cross(w=-1)   z(q=0)"]
for energy in ("lin", "quad"):
    for kappa in (0.5, 1, 2, 4, 16):
        r = g["rho_hist"](0.31, 0.68, kappa, energy, a[::-1])[::-1]; lr = np.log(r); la = np.log(a)
        w = -1 - np.gradient(lr, la)/3
        Or = g["W_R"]/0.68**2; H2 = 0.31*a**-3 + Or*a**-4 + r; q = (0.31*a**-3/2 + Or*a**-4 + (1 + 3*w)*r/2)/H2
        z = 1/a - 1; zc = z[np.where(np.diff(np.sign(w + 1)))[0]]; zq = z[np.where(np.diff(np.sign(q)))[0]]
        out.append(f" {energy:5s} {kappa:5.1f}  {w[-1]:+.3f}   {', '.join(f'{x:.2f}' for x in zc) or 'none':12s}   {zq[0]:.2f}")
out.append(" instant law (Claim 1): crossing exactly at q = 0")
txt = "\n".join(out); print(txt); open("iter1_wz.txt", "w").write(txt + "\n")
