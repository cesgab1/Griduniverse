"""
Is the strength of electromagnetism set by the grid sitting at a phase boundary (Multiple Point Principle)?
Step 1  geometry: continuum coupling of the grid's triangle action.  S = beta sum_t (1 - cos theta_t) ~ (beta/2) sum_t (F.a_t)^2
        -> e^2 = 1 / (beta * rho_A),  rho_A = sum_t A_t^2 / (6 V)   (hypercubic check: rho_A = 1, e^2 = 1/beta)
Step 2  phase boundary: compact U(1) Monte Carlo on the same 4D random grid, cold (ordered) and hot (random) starts,
        plaquette vs beta -> hysteresis / jump = the confinement-Coulomb transition beta_c
Step 3  critical fine-structure constant alpha_crit = kappa / (4 pi beta_c rho_A), kappa = measured coupling renormalisation
        (Wilson-loop fits near the transition), compared with the lattice-universal alpha_crit ~ 0.2
"""
import numpy as np, sys, os, json, time
os.environ["BOXL"] = "8"; start = sys.argv[1]; betas = [float(b) for b in sys.argv[2].split(",")]
sys.argv = ["x", "U1", "1.0", "1"]
src = open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/strong_on_grid.py").read().split("# ---- Wilson loops")[0]
exec(src)
DISP = np.array([l[2] for l in links])
u = TD[:, 0, None] * DISP[TL[:, 0]]; v = TD[:, 1, None] * DISP[TL[:, 1]]
A2 = 0.25 * ((u*u).sum(1) * (v*v).sum(1) - ((u*v).sum(1))**2)
rhoA = A2.sum() / (6 * L**4)
print(f"triangle geometry: mean area {np.sqrt(A2).mean():.3f}, rho_A = {rhoA:.4f}  ->  e^2 = {1/rhoA:.3f}/beta", flush=True)
out = []
for beta in (betas if start == "cold" else betas[::-1]):
    th = np.zeros(nl) if start == "cold" else rng.uniform(-np.pi, np.pi, nl)
    for s in range(int(os.environ.get("NTH", 60))): th = metro_u1(th, beta)
    pl = []
    for s in range(40): th = metro_u1(th, beta); pl.append(plaq_u1(th))
    out.append((beta, float(np.mean(pl)), float(np.std(pl) / np.sqrt(len(pl)))))
    print(f"{start} start  beta = {beta:.3f}: plaquette {np.mean(pl):.4f} ± {np.std(pl)/np.sqrt(len(pl)):.4f}   [{time.time()-t0:.0f}s]", flush=True)
json.dump(dict(rhoA=float(rhoA), scan=out), open(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/alpha_crit_{start}_{os.environ.get('TAG','a')}.json", "w"))
