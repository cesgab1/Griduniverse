"""Working model check: is the measured galaxy gap curve PREDICTIVE? Build it from half the galaxies (random split),
predict the other half's rotation from their visible matter only. Repeated over 200 random splits."""
import glob, os, numpy as np
KPC = 3.086e19
G = {}
for fn in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(fn, comments="#", ndmin=2); r, V, eV, Vg, Vd, Vb = d.T[:6]
    vb2 = Vg*np.abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2; ok = (V > 0) & (eV < 0.1*V) & (vb2 > 0)
    if ok.sum(): G[os.path.basename(fn)] = (np.log10(vb2[ok]*1e6/(r[ok]*KPC)), np.log10(V[ok]**2*1e6/(r[ok]*KPC)))
names = np.array(sorted(G)); rng = np.random.default_rng(7); res_pred, res_newton = [], []
for k in range(200):
    tr = rng.permutation(len(names))[:len(names)//2]; te = np.setdiff1d(np.arange(len(names)), tr)
    lg = np.concatenate([G[names[i]][0] for i in tr]); lo = np.concatenate([G[names[i]][1] for i in tr])
    edges = np.arange(-12.6, -8.0, 0.2); c, m = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (lg >= a) & (lg < b)
        if s.sum() >= 10: c.append((a+b)/2); m.append(np.median(lo[s] - lg[s]))
    lgt = np.concatenate([G[names[i]][0] for i in te]); lot = np.concatenate([G[names[i]][1] for i in te])
    res_pred.append(np.std(lot - (lgt + np.interp(lgt, c, m)))); res_newton.append(np.median(np.abs(lot - lgt)))
out = [f"Hold-out test, 200 random half/half splits of {len(names)} galaxies (curve built on one half, applied to the other):",
       f"  prediction error with the measured gap curve: {np.mean(res_pred):.3f} dex (x{10**np.mean(res_pred):.2f}) typical",
       f"  Newton with visible matter alone: off by {np.mean(res_newton):.3f} dex (x{10**np.mean(res_newton):.2f}) typical (median)"]
txt = "\n".join(out); print(txt); open("holdout.txt", "w").write(txt + "\n")
