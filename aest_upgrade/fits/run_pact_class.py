"""Upgraded AeST in CLASS (our patched classy): P-ACT (Planck low-l + Planck cut + ACT DR6) + DESI DR2 BAO + DES-Dovekie SN.
Models (env AEST_BETA, AEST_EX): AeST(Cosh, paper values) + Λ; + grid-tension dark energy β; + energy exchange with the condensate."""
import sys, os, json, time
sys.path.insert(0, "/home/claude/act_lite")
from cobaya.run import run
tag = sys.argv[1]; start = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
beta = float(os.environ.get("AEST_BETA", 0)); ex = os.environ.get("AEST_EX", "no")
extra = {"gauge": "newtonian", "non_linear": "halofit", "l_max_scalars": 7000, "N_ur": 2.0328, "N_ncdm": 1, "m_ncdm": 0.06,
         "omega_cdm": 0.0, "fluid_equation_of_state": "AEST", "use_ppf": "no",
         "aest_KB": 0.5, "aest_K2": 7.5e3, "aest_Q0": 0.1, "aest_Z0": 1e-9, "aest_de_beta": beta, "aest_exchange": ex}
def P(name, lo, hi, ref, prop): return {"prior": {"min": lo, "max": hi}, "ref": start.get(name, ref), "proposal": prop}
params = {"logA": dict(P("logA", 2.9, 3.2, 3.05, 0.004), drop=True), "A_s": {"value": "lambda logA: 1e-10*np.exp(logA)"},
          "n_s": P("ns", 0.9, 1.02, 0.968, 0.003), "H0": P("H0", 55, 85, 68.0, 0.3), "omega_b": P("ombh2", 0.02, 0.025, 0.0224, 0.0001),
          "omch2": dict(P("omch2", 0.08, 0.16, 0.119, 0.0008), drop=True),
          "Omega_fld": {"value": "lambda omch2, H0: omch2/(H0/100.)**2"},
          "tau_reio": P("tau", 0.01, 0.12, 0.058, 0.004),
          "A_act": {"prior": {"min": 0.5, "max": 1.5}, "ref": start.get("A_act", 1.0), "proposal": 0.002},
          "P_act": {"prior": {"min": 0.9, "max": 1.1}, "ref": start.get("P_act", 1.0), "proposal": 0.003}}
info = {"packages_path": "/home/claude/cobaya_packages",
        "likelihood": {"planck_2018_lowl.TT": None, "planck_2018_lowl.EE_sroll2": None,
                       "act_dr6_cmbonly.PlanckActCut": {"dataset_params": {"use_cl": "tt te ee", "lmin_cuts": "0 0 0", "lmax_cuts": "1000 600 600"},
                                                        "params": {"A_planck": {"value": "lambda A_act: A_act"}}},
                       "act_dr6_cmbonly": {"lmax_theory": 7000}, "bao.desi_dr2": None, "sn.desdovekie": None},
        "prior": {"cal_dip_prior": "lambda A_act: stats.norm.logpdf(A_act, loc = 1.0, scale = 0.003)"},
        "theory": {"classy": {"extra_args": extra}}, "params": params,
        "sampler": {"minimize": {"method": "bobyqa", "best_of": 1, "max_evals": 300}} if "EVAL" not in os.environ else {"evaluate": None},
        "output": f"/home/claude/fullfit/pact_class/{tag}", "force": True}
os.makedirs("/home/claude/fullfit/pact_class", exist_ok=True)
t0 = time.time(); upd, sampler = run(info)
f = f"/home/claude/fullfit/pact_class/{tag}" + (".minimum.txt" if "EVAL" not in os.environ else ".1.txt")
cols = open(f).readline().lstrip("#").split(); vals = open(f).readlines()[1].split(); m = dict(zip(cols, map(float, vals)))
json.dump({"chi2": m.get("chi2"), "parts": {k: v for k, v in m.items() if k.startswith("chi2__")},
           "params": {k: m[k] for k in ("H0","n_s","omega_b","omch2","tau_reio","A_act","logA") if k in m}, "minutes": (time.time()-t0)/60},
          open(f"/home/claude/fullfit/pact_class/{tag}.json", "w"), indent=1)
print("done", (time.time()-t0)/60, "min")
