import numpy as np
exec(open("/tmp/claude-0/-home-claude-griduniverse/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/kb_direct.py").read().split("ref = None")[0])
def run(kb, prec):
    c = Class(); c.set({**base, "aest_KB": kb, **prec}); c.compute(); tt = c.lensed_cl(2500)["tt"][2:2501]; c.struct_cleanup(); c.empty(); return tt
P0 = {}; P1 = {"tol_perturbations_integration": 1e-7, "perturbations_sampling_stepsize": 0.02, "start_small_k_at_tau_c_over_tau_h": 0.0004, "start_large_k_at_tau_h_over_tau_k": 0.02}
ref = run(0.5, P1)
for kb in (0.1, 0.07):
    a = run(kb, P0)/ref-1; b = run(kb, P1)/ref-1
    print(f"K_B={kb}: default precision l=1000 {a[998]*100:+.3g}%  max {np.max(abs(a))*100:.3g}% | high precision l=1000 {b[998]*100:+.3g}%  max {np.max(abs(b))*100:.3g}%", flush=True)
