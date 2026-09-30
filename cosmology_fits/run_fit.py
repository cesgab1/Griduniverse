import sys, os, json, time
sys.path.insert(0, "/home/claude/fullfit")
from cobaya.run import run
model, sn = sys.argv[1], sys.argv[2]
refine = len(sys.argv) > 3 and sys.argv[3] in ("refine", "fromlaw", "mnu", "gboost", "gfluid")
g_mode = len(sys.argv) > 3 and sys.argv[3] in ("gboost", "gfluid"); gf_mode = g_mode and sys.argv[3] == "gfluid"
GD = os.environ.get("GDELTA", "0") + ("f" + os.environ.get("GCS2", "1") if gf_mode else "") + ("_at" + os.environ["GAT"] if "GAT" in os.environ else "") + ("_" + os.environ["GFORM"] if "GFORM" in os.environ else "")
MNU = float(os.environ.get("MNU", 0.06)); mnu_mode = len(sys.argv) > 3 and sys.argv[3] == "mnu"
start_tag = "_fromlaw" if (len(sys.argv) > 3 and sys.argv[3] == "fromlaw") else ""
import pandas as pd
def read_min(path):
    cols = open(path).readline().lstrip("#").split()
    vals = open(path).readlines()[1].split()
    return dict(zip(cols, map(float, vals)))
extra = {"lens_potential_accuracy": 1, "num_massive_neutrinos": 1, "nnu": 3.044, "halofit_version": "mead2020"}
params = {"logA": {"prior": {"min": 2.9, "max": 3.2}, "ref": 3.045, "proposal": 0.01, "drop": True},
          "As": {"value": "lambda logA: 1e-10*np.exp(logA)"},
          "ns": {"prior": {"min": 0.9, "max": 1.02}, "ref": 0.968, "proposal": 0.004},
          "H0": {"prior": {"min": 55, "max": 85}, "ref": 67.5, "proposal": 0.5},
          "ombh2": {"prior": {"min": 0.02, "max": 0.025}, "ref": 0.0222, "proposal": 0.00015},
          "omch2": {"prior": {"min": 0.08, "max": 0.16}, "ref": 0.119, "proposal": 0.001},
          "tau": {"prior": {"min": 0.01, "max": 0.12}, "ref": 0.058, "proposal": 0.006},
          "mnu": MNU}
vc_mode = "VC_ME" in os.environ
if model == "LCDM":
    theory = {("gboost_fluid.CAMBGF" if gf_mode else "gboost_camb.CAMBG") if g_mode else ("vc_camb.CAMBVC" if vc_mode else "camb"): {"extra_args": extra}}
elif model == "w0wa":
    theory = {"camb": {"extra_args": dict(extra, dark_energy_model="ppf")}}
    params["w"] = {"prior": {"min": -3, "max": 1}, "ref": -0.8, "proposal": 0.05}
    params["wa"] = {"prior": {"min": -3, "max": 2}, "ref": -0.6, "proposal": 0.2}
else:   # LAW with beta
    theory = {(("gboost_fluid.CAMBLawGF" if gf_mode else "gboost_camb.CAMBLawG") if g_mode else "law_camb.CAMBLaw"): ({"extra_args": extra} if gf_mode else {"beta": float(model.split("_")[1]), "extra_args": extra})}
info = {"packages_path": "/home/claude/cobaya_packages",
        "likelihood": {"planck_2018_lowl.TT": None, "planck_2018_lowl.EE_sroll2": None,
                       "planck_NPIPE_highl_CamSpec.TTTEEE": None, "bao.desi_dr2": None, sn: None},
        "theory": theory, "params": params,
        "sampler": {"minimize": {"method": "bobyqa", "best_of": 1 if refine else 2, "max_evals": "3000d"}},
        "output": f"/home/claude/fullfit/out/{model}_{sn}" + ((f"_mnu{MNU}" if mnu_mode else (f"_g{GD}" if g_mode else (start_tag or "_refine"))) if refine else "") + ("_me" + os.environ["VC_ME"] if "VC_ME" in os.environ else ""), "force": True}
if refine:
    start = read_min(f"/home/claude/fullfit/out/{model}_{sn}" + ("_refine" if (mnu_mode or g_mode) else ("_lawstart" if start_tag else "")) + ".minimum.txt")
    for k, v in info["params"].items():
        if isinstance(v, dict) and "prior" in v and k in start:
            v["ref"] = start[k]; v["proposal"] = v.get("proposal", 0.01) / 3
    nuis = {k: start[k] for k in ["A_planck","amp_143","amp_217","amp_143x217","n_143","n_217","n_143x217","calTE","calEE"]}
    from cobaya.component import get_component_class
    dflt = get_component_class("planck_NPIPE_highl_CamSpec.TTTEEE", kind="likelihood").get_defaults()["params"]
    for k, v in nuis.items():
        d = dict(dflt[k]); d["ref"] = v
        if "proposal" in d: d["proposal"] = d["proposal"] / 3
        info["params"][k] = d
t0 = time.time()
upd, sampler = run(info)
tag = ((f"_mnu{MNU}" if mnu_mode else (f"_g{GD}" if g_mode else (start_tag or "_refine"))) if refine else "") + ("_me" + os.environ["VC_ME"] if vc_mode else "")
m = read_min(f"/home/claude/fullfit/out/{model}_{sn}{tag}.minimum.txt")
res = {"model": model, "sn": sn, "refine": refine, "minutes": round((time.time() - t0) / 60, 1), "chi2": m["chi2"],
       "parts": {k: v for k, v in m.items() if k.startswith("chi2__") and k.count("__") == 1 and "." in k}}
res["params"] = {k: m[k] for k in ["ns", "H0", "ombh2", "omch2", "tau", "w", "wa"] if k in m}
json.dump(res, open(f"/home/claude/fullfit/out/{model}_{sn}{tag}.json", "w"), indent=1)
print(json.dumps(res, indent=1))
