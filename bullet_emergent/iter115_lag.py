"""Iteration 115: delayed emergent-gravity response in the Bullet Cluster (PREREG_115.md)."""
import os, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "iter114_bullet_eg.py")).read().split("out = [")[0]; exec(src)
KPC_PER = 1.0227e-3                                     # kpc per (km/s * Myr)
Mm = gas[0][1] + st[0][1]; Ms = gas[1][1] + st[1][1]
e = (P["sub BCG"] - P["main BCG"]); e = e/np.linalg.norm(e)
vrel, T = 2700.0, 150.0
vs, vm = vrel*Mm/(Mm + Ms), vrel*Ms/(Mm + Ms)
def lagged(tau, f):
    out = []
    for (c_, M, a), d, v in [(gas[0][0], gas[0][1], gas[0][2]), (gas[1][0], gas[1][1], gas[1][2])] and []: pass
    for i, (sign, v) in enumerate([(-1, vm), (+1, vs)]):
        cs, Ms_, as_ = st[i]; cg, Mg, ag = gas[i]
        out.append((cs - sign*e*v*tau*KPC_PER, Ms_, as_))
        back = f*v*min(tau, T) + v*max(tau - T, 0)
        out.append((cg - sign*e*back*KPC_PER, Mg, ag))
    return out
h = 20.0; gx = np.arange(-900, 1401, h); gy = np.arange(-600, 601, h); gz = np.arange(-3000, 3001, h)
res = [f"Iteration 115 -- lagged emergent pull; v_sub {vs:.0f}, v_main {vm:.0f} km/s along the axis; T = {T:.0f} Myr", ""]
saved = list(blobs)
for f in (0.5, 0.0, 1.0):
    res.append(f"gas speed after passage = {f} x galaxies'")
    best = None
    for tau in range(0, 401, 25):
        blobs[:] = lagged(tau, f)
        X, Y, Z, g, rho = fields(gx, gy, gz); gm = np.sqrt(g[0]**2 + g[1]**2 + g[2]**2); fac = np.sqrt(aV/gm)
        div = sum(np.gradient(fac*g[i], h, axis=i) for i in range(3)); Sig = (-div/(4*np.pi*Gk)).sum(axis=2)*h
        blobs[:] = saved                                   # Newton + tip from matter where it is NOW
        import io, contextlib
        with contextlib.redirect_stdout(io.StringIO()): lines, k = score("", gx, gy, gz, h, Sig_ph=Sig)
        m, s = k["main BCG"] - k["main plasma"], k["sub BCG"] - k["sub plasma"]
        zc = np.hypot((dM - m)/emain, (dS - s)/esub)
        res.append(f"   tau {tau:3d} Myr: main {m:+.3f} sub {s:+.3f} -> combined {zc:.1f} sigma; {lines[3].strip()}")
        if best is None or zc < best[0]: best = (zc, tau)
    res.append(f"   BEST: tau = {best[1]} Myr, {best[0]:.1f} sigma"); res.append("")
txt = "\n".join(res); print(txt); open(os.path.join(here, "iter115_lag.txt"), "w").write(txt + "\n")
