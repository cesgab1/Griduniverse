"""
ITERATION 35 (QG step 2): can stacked layers shift the grid's smallest-scale scaling from z = 3 to z = 3 - eta?
Estimate (written before computing): N layers ('species') saturating the species bound (N G Lambda_cut^2 ~ 1) make the one-loop
correction to the scaling exponent of order eta ~ c / (16 pi^2), with c an O(1) coefficient we CANNOT compute here (it needs
the full Horava loop calculation). Its SIGN is also not determined at this level.
Criterion: order-of-magnitude consistency only -- does c in 0.5-3 bracket the needed eta = 0.013 (ACT, n_s 0.974) to 0.017
(Planck, n_s 0.965)? No value of c is picked.
"""
import numpy as np
base = 1/(16*np.pi**2)
out = ["ITERATION 35: watermark tilt from stacked layers (order-of-magnitude only)", "",
       f"natural size eta ~ c/(16 pi^2) = c x {base:.4f}", " c      eta       n_s if the shift lowers z (red)    n_s if it raises z (blue)"]
for c in (0.5, 1, 2, 3):
    eta = c*base; ns_red = 1 - 2*eta/(1 - eta); ns_blue = 1 + 2*eta/(1 + eta)
    out.append(f" {c:3.1f}   {eta:.4f}        {ns_red:.4f}                            {ns_blue:.4f}")
out += ["", "needed: eta = 0.013 (ACT) - 0.017 (Planck), red direction (z just below 3).",
        "Reading: the natural one-loop size from layers saturating the species bound (0.003-0.019) brackets what is needed, so the",
        "idea is not off by orders of magnitude. But the coefficient AND the sign are uncomputed here; matching n_s would mean",
        "choosing c ~ 1.7-2.2 and the red sign by hand -> no prediction yet (accommodation). Needs the real loop calculation."]
txt = "\n".join(out); print(txt); open("iter35_watermark_eta.txt", "w").write(txt + "\n")
