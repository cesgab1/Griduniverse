"""GR track G0c (Coalesce): what Special Relativity says about the two places GR stops, and how they connect. Numbers only."""
import numpy as np
c, H0 = 299792.458, 67.4
Rh = c/H0/1000                                                      # Hubble radius, Gpc
Vvis = 4/3*np.pi*14.1**3
out = ["GR track G0c: Special Relativity and the two gaps", ""]
for ok in (0.0012, 0.0023, 0.0034):
    R = Rh/np.sqrt(ok); Vmin = 0.9427*R**3                          # smallest finite hyperbolic space (Weeks manifold), volume 0.9427 R^3
    out.append(f"Omega_k = {ok:.4f}: curvature radius {R:.0f} Gpc; smallest FINITE open-curved space has volume {Vmin:.1e} Gpc^3 = {Vmin/Vvis:.0f} x visible;"
               f" its shortest loop ~{0.58*R:.0f} Gpc (matched circles need < 27.7 Gpc)")
txt = "\n".join(out); print(txt); open("G0c_special_relativity.txt", "w").write(txt + "\n")
