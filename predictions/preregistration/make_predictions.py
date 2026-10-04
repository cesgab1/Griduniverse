"""
Pre-registered predictions, Grid Universe (frozen 2026-10-03, before DESI DR3 / Euclid / Simons Observatory results).
Computes, for the best fits to DESI DR2 BAO + Planck 2018 distance priors + each supernova set:
  w(z) curves, w = -1 crossing redshift, bin-averaged w, H(z)/H_LCDM, and BAO distances D_M/r_d, D_H/r_d at DESI tracer
  redshifts, for: LCDM, CLAIM 1 (instant Tension-Rate Law, beta = 1/2), TOY (memory kappa = 3 from Hubble friction,
  quadratic energy, VCDM reading). No parameter is tuned here: all come from the stored fits.
"""
import numpy as np, json, os, re
os.environ["ONLY"] = "NONE"
FIT = "/home/claude/griduniverse/cosmology_fits/fit_law.py"
TOY = "/home/claude/griduniverse/quantum_gravity/toy_grid/iter1_memory.py"
res = {}; zs_bins = [(0.1, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.1), (1.1, 1.6), (1.6, 2.1)]
z_bao = [0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330]
for SN in ("PANTHEON", "DESY5", "UNION3"):
    os.environ["SNSET"] = SN
    g = {"__name__": "x", "__file__": TOY}; exec(open(TOY).read().split("# ---------- Part A")[0], g); fl = g["g"]
    txt = open(f"/home/claude/griduniverse/predictions/fit_{SN}_beta0.5_rsfix.txt").read()
    par = lambda m: [float(x) for x in re.search(m + r"\s+chi2.*params \[([^\]]+)\]", txt).group(1).split()]
    toy = json.load(open(f"/home/claude/griduniverse/quantum_gravity/toy_grid/iter1_fits_{SN}_k3.json"))["quad kappa=3.0"][1]
    E_orig = fl["E_of_z"]; E_toy = g["make_E"](3.0, "quad"); ZG = fl["ZG"]
    out = {}
    for name, p, Ef, model in (("LCDM", par("LCDM"), E_orig, "LCDM"), ("CLAIM1", par("LAW"), E_orig, "LAW"), ("TOY", toy, E_toy, "TOY")):
        fl["E_of_z"] = Ef; Om, h, wb = p; E = Ef(model, Om, h, [])
        Or = fl["W_R"]/h**2; rde = E**2 - Om*(1 + ZG)**3 - Or*(1 + ZG)**4
        m = ZG <= 3; z = ZG[m]; lr = np.log(np.clip(rde[m], 1e-30, None))
        w = -1 + np.gradient(lr, np.log(1 + z))/3 if name != "LCDM" else -np.ones_like(z)
        cross = z[np.where(np.diff(np.sign(w + 1)))[0]] if name != "LCDM" else []
        o = fl["observables"](model, p)
        bao = {f"{zz}": [float(o["DM"](zz)/o["rd"]), float(o["DH"](zz)/o["rd"])] for zz in z_bao}
        wb_ = {f"{a}-{b}": float(np.mean(w[(z >= a) & (z < b)])) for a, b in zs_bins}
        out[name] = dict(params=p, w0=float(w[0]), crossing=[float(c) for c in cross], w_bins=wb_, bao=bao,
                         H0=100*h, Om=Om, z=z[::25].tolist(), w=w[::25].tolist(), E=E[m][::25].tolist())
    fl["E_of_z"] = E_orig
    for name in ("CLAIM1", "TOY"):
        out[name]["H_over_LCDM"] = (np.array(out[name]["E"])*out[name]["H0"]/(np.array(out["LCDM"]["E"])*out["LCDM"]["H0"])).tolist()
    res[SN] = out
json.dump(res, open("predictions.json", "w"), indent=1)
lines = []
for SN in res:
    for name in ("CLAIM1", "TOY", "LCDM"):
        r = res[SN][name]
        lines.append(f"{SN:9s} {name:6s} H0 {r['H0']:.2f} Om {r['Om']:.4f} w0 {r['w0']:+.3f} crossing {r['crossing']}  bins " +
                     " ".join(f"{k}:{v:+.3f}" for k, v in r["w_bins"].items()))
open("predictions_table.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
