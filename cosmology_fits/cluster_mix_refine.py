exec(open("cluster_mix.py").read().split("out = []")[0])
import warnings; warnings.filterwarnings("ignore")
fs = np.round(np.arange(0.80, 1.001, 0.025), 3); res = []
for f in fs:
    rr = []
    for M, b, s in zip(Ms, base, sig):
        d = thr(M*Msun, f, True); dc = 1.686*d/b; rr.append(erfc(dc/np.sqrt(2)/s)/erfc(1.686/np.sqrt(2)/s))
    chi, med = pull_chi(f); res.append((f, rr, med, chi))
    print(f"f {f:.3f}: counts " + " ".join(f"×{x:6.3f}" for x in rr) + f" | pull ×{med:.2f}")
# where each requirement is met
fcount = [np.interp(0, [np.log(x[1][k]) for x in res], fs) for k in range(3)]
print("counts match (×1) at f =", ", ".join(f"{v:.3f}" for v in fcount), "(M = 1e14, 3e14, 1e15)")
ff = np.linspace(0, 1, 1001); meds = [pull_chi(x)[1] for x in ff]; chis = [pull_chi(x)[0] for x in ff]
print(f"per-cluster pull: median matches at f = {np.interp(1, meds, ff):.3f}; best shape fit at f = {ff[np.argmin(chis)]:.3f}")
f9 = 0.9; gN = gbar*(1+f9*Ox/Ob); gm = gN*nu(gN/a0)
print("f=0.9 model/measured by g_bar:", " ".join(f"{x:.0e}:{y:.2f}" for x, y in zip(gbar[::3], (gm/gmeas)[::3])))
# eROSITA: S8 = 0.86 ± 0.01 -> tolerance on counts: abundance change for ±1σ in σ8 at 3e14
s = sig[1]; nu0 = 1.686/s
for ds in (-0.012, 0.012):
    print(f"counts tolerance from eROSITA ±1σ (σ8 ±1.2%): ×{erfc(nu0/(1+ds)/np.sqrt(2))/erfc(nu0/np.sqrt(2)):.2f} at 3e14")
