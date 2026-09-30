"""P-ACT combination (ACT DR6 recommended): Planck low-l + Planck plik-lite cut at l<1000 TT / <600 TE,EE + ACT DR6 CMB-only
+ DESI DR2 BAO + DES-Dovekie SN. Model via gboost_fluid (GDELTA=0 -> LCDM; GFORM recoil/waves; GAT; GCS2)."""
import sys, os, json, time
sys.path.insert(0, "/home/claude/fullfit"); sys.path.insert(0, "/home/claude/act_lite")
from cobaya.run import run
tag = sys.argv[1]; start = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
extra = {"lens_potential_accuracy": 4, "lAccuracyBoost": 1.1, "kmax": 10, "nonlinear": True,
         "num_massive_neutrinos": 1, "nnu": 3.044, "halofit_version": "mead2020"}
def P(name, lo, hi, ref, prop): return {"prior": {"min": lo, "max": hi}, "ref": start.get(name, ref), "proposal": prop}
params = {"logA": dict(P("logA", 2.9, 3.2, 3.05, 0.004), drop=True), "As": {"value": "lambda logA: 1e-10*np.exp(logA)"},
          "ns": P("ns", 0.9, 1.02, 0.968, 0.003), "H0": P("H0", 55, 85, 68.0, 0.3), "ombh2": P("ombh2", 0.02, 0.025, 0.0224, 0.0001),
          "omch2": P("omch2", 0.08, 0.16, 0.119, 0.0008), "tau": P("tau", 0.01, 0.12, 0.058, 0.004), "mnu": 0.06,
          "A_act": {"prior": {"min": 0.5, "max": 1.5}, "ref": start.get("A_act", 1.0), "proposal": 0.002}}
if os.environ.get("OMK") == "free": params["omk"] = {"prior": {"min": -0.05, "max": 0.05}, "ref": 0.0, "proposal": 0.001}
info = {"packages_path": "/home/claude/cobaya_packages",
        "likelihood": {"planck_2018_lowl.TT": None, "planck_2018_lowl.EE_sroll2": None,
                       "act_dr6_cmbonly.PlanckActCut": {"dataset_params": {"use_cl": "tt te ee", "lmin_cuts": "0 0 0", "lmax_cuts": "1000 600 600"},
                                                        "params": {"A_planck": {"value": "lambda A_act: A_act"}}},
                       "act_dr6_cmbonly": None, "bao.desi_dr2": None, "sn.desdovekie": None},
        "prior": {"cal_dip_prior": "lambda A_act: stats.norm.logpdf(A_act, loc = 1.0, scale = 0.003)"},
        "theory": {("vc_camb.CAMBVC" if "VC_ME" in os.environ or "VC_ALPHA" in os.environ else "gboost_fluid.CAMBGF"): {"extra_args": extra}}, "params": params,
        "sampler": {"minimize": {"method": "bobyqa", "best_of": 1, "max_evals": "3000d"}} if "EVAL" not in os.environ else {"evaluate": None},
        "output": f"/home/claude/fullfit/pact/{tag}", "force": True}
t0 = time.time(); upd, sampler = run(info)
if "EVAL" not in os.environ:
    f = f"/home/claude/fullfit/pact/{tag}.minimum.txt"; cols = open(f).readline().lstrip("#").split(); vals = open(f).readlines()[1].split()
    m = dict(zip(cols, map(float, vals)))
    json.dump({"chi2": m["chi2"], "parts": {k: v for k, v in m.items() if k.startswith("chi2__")}, "params": {k: m[k] for k in ("H0","ns","ombh2","omch2","tau","A_act") if k in m},
               "minutes": (time.time()-t0)/60}, open(f"/home/claude/fullfit/pact/{tag}.json", "w"), indent=1)
print("done", (time.time()-t0)/60, "min")
