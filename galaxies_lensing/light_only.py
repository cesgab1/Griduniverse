"""
Galaxy spin from light alone: blind prediction of rotation curves, no per-galaxy tuning.
Inputs per galaxy: only the stars (3.6-micron light x fixed mass-to-light 0.5 disk / 0.7 bulge, from stellar-population models)
and the measured gas. One universal number a0 -- and to keep it honest, a0 is fitted on HALF the galaxies and
the predictions are made for the OTHER half (then swapped). Compared with the measured speeds, point by point.
Benchmarks: Newton with visible matter only (no slack, no dark matter).
"""
import numpy as np, glob, os, json
from scipy.optimize import minimize_scalar
kpc = 3.0857e19
gals = []
for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat")):
    d = np.loadtxt(f, comments="#"); d = d[None] if d.ndim == 1 else d
    r, V, eV, Vg, Vd, Vb = d.T[:6]
    ok = (V > 0) & (r > 0)
    if ok.sum() < 5: continue
    gals.append(dict(name=os.path.basename(f)[:-11], r=r[ok], V=V[ok], eV=eV[ok], Vg=Vg[ok], Vd=Vd[ok], Vb=Vb[ok]))
def vpred(g, a0, slack=True):
    gN = (g["Vg"]*np.abs(g["Vg"]) + 0.5*g["Vd"]**2 + 0.7*g["Vb"]**2)*1e6/(g["r"]*kpc)
    gN = np.clip(gN, 1e-16, None)
    nu = 1/(1 - np.exp(-np.sqrt(gN/a0))) if slack else 1.0
    return np.sqrt(gN*nu*g["r"]*kpc)/1e3
def score(G, a0):
    return np.concatenate([np.log10(vpred(g, a0)/g["V"]) for g in G])
rng = np.random.default_rng(1); idx = rng.permutation(len(gals)); A, B = [gals[i] for i in idx[::2]], [gals[i] for i in idx[1::2]]
res = []
for train, test in ((A, B), (B, A)):
    a0 = 10**minimize_scalar(lambda x: np.sum(score(train, 10**x)**2), bounds=(-11, -9), method="bounded").x
    for g in test:
        vp = vpred(g, a0); vn = vpred(g, a0, slack=False)
        res.append(dict(name=g["name"], a0=a0, rms=float(np.sqrt(np.mean((vp/g["V"]-1)**2))), outer=float(vp[-1]/g["V"][-1]-1),
                        newton_outer=float(vn[-1]/g["V"][-1]-1), vmax=float(g["V"].max()),
                        within_err_frac=float(np.mean(np.abs(vp-g["V"]) < np.maximum(g["eV"], 0.05*g["V"])))))
    print(f"a0 fitted on {len(train)} galaxies: {a0:.3e} m/s^2  -> predicting the other {len(test)}")
rms = np.array([x["rms"] for x in res]); outer = np.array([x["outer"] for x in res]); nout = np.array([x["newton_outer"] for x in res])
print(f"\nBlind predictions for {len(res)} galaxies (speed at every measured radius, from light + gas only):")
print(f"  typical point-by-point error: median {np.median(rms)*100:.0f}%   (16-84%: {np.percentile(rms,16)*100:.0f}%-{np.percentile(rms,84)*100:.0f}%)")
print(f"  outermost measured speed: median error {np.median(outer)*100:+.1f}%, within 10% for {np.mean(abs(outer)<0.10)*100:.0f}%, within 20% for {np.mean(abs(outer)<0.20)*100:.0f}% of galaxies")
print(f"  Newton with visible matter only, outermost speed: median error {np.median(nout)*100:+.0f}% (within 20% for {np.mean(abs(nout)<0.2)*100:.0f}%)")
print(f"  share of individual points matching within their error bar (or 5%): {np.mean([x['within_err_frac'] for x in res])*100:.0f}%")
for lo, hi, lab in ((0, 60, "dwarfs (<60 km/s)"), (60, 150, "medium"), (150, 400, "giants (>150 km/s)")):
    m = np.array([lo <= x["vmax"] < hi for x in res])
    print(f"  {lab:20s} n={m.sum():3d}: median rms error {np.median(rms[m])*100:.0f}%, outer speed median {np.median(outer[m])*100:+.1f}%")
bad = sorted(res, key=lambda x: -x["rms"])[:6]
print("  worst:", ", ".join(f"{x['name']} ({x['rms']*100:.0f}%)" for x in bad))
json.dump(res, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/light_only.json", "w"), indent=1)
