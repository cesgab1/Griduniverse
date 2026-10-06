"""Iteration 113 POST-HOC (not pre-registered): galaxy rule (MOND, a0 = 1.2e-10) on the same Coma profiles, for comparison."""
import os, numpy as np
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "iter113_emergent.py")).read().split("radii = [")[0]; exec(src)
a0 = 1.2e-10; out = ["POST-HOC: MOND (simple nu, a0 1.2e-10) vs emergent gravity, same visible matter; ratios to hydrostatic mass"]
for kind in ("gas-shape", "hernquist"):
    MB = Mgas + stars(kind); gN = G*MB/r**2; nu = 1/(1 - np.exp(-np.sqrt(gN/a0)))
    MD = np.sqrt(np.clip(c*H0*r**2/(6*G)*np.gradient(MB*r, r), 0, None))
    out.append(f" stars {kind}: " + ", ".join(f"{R} kpc MOND {np.interp(R*KPC, r, MB*nu/Mmeas):.2f} / EG {np.interp(R*KPC, r, (MB+MD)/Mmeas):.2f}" for R in [100, 300, 500, 1000, 1500, 2500]))
Mvis = np.interp(2500*KPC, r, Mgas + stars("gas-shape"))/MS
out.append(f" verified visible mass within 2.5 Mpc {Mvis:.2e} Msun = {Mvis/1.69e15:.3f} of hydrostatic -> short x{1.69e15/Mvis:.1f} (iteration 112 used 3e14 -> x5.6)")
txt = "\n".join(out); print(txt); open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "iter113_posthoc_mond.txt"), "w").write(txt + "\n")
