"""H0 with NO CMB light at all: DESI DR2 BAO + BBN baryon density (from deuterium, Cooke et al.: omega_b = 0.02218 +- 0.00055)
(+ optionally DES-Dovekie supernova shape). Sound-horizon ruler computed by CAMB from ordinary early physics."""
import sys; sys.path.insert(0, "/home/claude/fullfit")
from cobaya.run import run
for lik in (["bao.desi_dr2"], ["bao.desi_dr2", "sn.desdovekie"]):
    info = {"packages_path": "/home/claude/cobaya_packages",
            "likelihood": {l: None for l in lik},
            "theory": {"camb": {"extra_args": {"num_massive_neutrinos": 1, "nnu": 3.044}}},
            "params": {"H0": {"prior": {"min": 50, "max": 90}, "ref": 68, "proposal": 0.5},
                       "ombh2": {"prior": {"dist": "norm", "loc": 0.02218, "scale": 0.00055}, "ref": 0.0222, "proposal": 0.0003},
                       "omch2": {"prior": {"min": 0.05, "max": 0.2}, "ref": 0.12, "proposal": 0.003},
                       "mnu": 0.06, "As": 2.1e-9, "ns": 0.965, "tau": 0.055},
            "sampler": {"mcmc": {"Rminus1_stop": 0.02, "max_tries": 10000, "learn_proposal": True}},
            "output": f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/noCMB_{len(lik)}", "force": True}
    upd, s = run(info)
    import numpy as np
    prod = s.products()["sample"]; h = prod["H0"].values; w = prod["weight"].values
    m = np.average(h, weights=w); sd = np.sqrt(np.average((h-m)**2, weights=w))
    print(f"\n### {' + '.join(lik)} + BBN (no CMB): H0 = {m:.2f} +- {sd:.2f}   (SH0ES 73.04 +- 1.04: {abs(73.04-m)/np.hypot(sd,1.04):.1f} sigma apart)\n")
