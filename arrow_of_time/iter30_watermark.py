"""
ITERATION 30 (pre-registration part): where did the 'watermark' (primordial ripples) come from in the grid?
Measured targets (values looked up AFTER committing this file): tilt n_s (Planck 2018: 0.965 +/- 0.004), running alpha_s,
amplitude A_s = 2.1e-9.
Two grid-native candidates, no others tried:
 W1 cell jitter: independent random nudges per cell / thermal noise of the jostling radiation, set locally after the start.
    Local, uncorrelated noise has no correlations beyond the horizon -> white-noise density field -> P_zeta(k) ~ k^3,
    i.e. n_s = 4 (strongly 'blue'). Prediction: EXCLUDED (Planck: n_s ~ 0.965).
 W2 Horava-type scaling at the cell scale (BORROWED: Mukohyama, JCAP 0906:001, arXiv:0904.2190). Our base has a preferred time
    slicing; at the shortest scales the dispersion is omega = p^z / M^(z-1) with z = 3 (the power-counting value in 3 space
    dimensions). A field freezes when omega = H; its spectrum is P ~ p_f^3 / H with p_f = (H M^(z-1))^(1/z):
       P ~ H^(3/z - 1) M^(3(z-1)/z).   For z = 3 the H-dependence cancels: EXACTLY scale invariant (n_s = 1) in ANY background.
    For z = 3 - eta, in the radiation era: n_s - 1 = -2 (3 - z)/(z - 2)  (derived below and checked numerically).
    Predictions: (a) pure z = 3 gives n_s = 1 -> EXCLUDED if data say n_s < 1 at > 5 sigma.
                 (b) matching the measured tilt needs eta = 3 - z ~ 0.017 (one fitted number = accommodation);
                     it then PREDICTS no running: |alpha_s| < 1e-3 (constant z in one era).
                 (c) amplitude: delta phi ~ M/(2 pi), but conversion to curvature needs a curvaton-type step with a free
                     efficiency -> amplitude NOT predicted.
"""
import numpy as np

def tilt_numeric(z, era_w=1/3, M=1.0):
    """freeze-out numerically: background H ~ a^(-3(1+w)/2); for a set of comoving k find a_freeze where k/a = (H M^(z-1))^(1/z),
    P(k) = p_f^3 / H at freeze; return d ln P / d ln k"""
    a = np.geomspace(1e-30, 1e-10, 200000)
    H = a**(-1.5*(1 + era_w))*1e-50
    pf = (H*M**(z - 1))**(1/z)
    kk = a*pf                                    # comoving k freezing at a
    P = pf**3/H
    sel = slice(50000, 150000)
    return np.polyfit(np.log(kk[sel]), np.log(P[sel]), 1)[0]

out = ["ITERATION 30: the grid's watermark -- tilt of the primordial ripples (criteria committed before the data lookup)", "",
       "W1 cell jitter (local, uncorrelated): P ~ k^3 -> n_s = 4 (no free choice).", "",
       "W2 Horava scaling, freeze-out tilt (numeric vs formula), radiation era:"]
for z in (3.0, 2.99, 2.983, 2.97, 3.01):
    num = tilt_numeric(z); formula = -2*(3 - z)/(z - 2)
    out.append(f"   z = {z:6.3f}: n_s - 1 = {num:+.4f} (numeric)   {formula:+.4f} (formula)")
out += ["   matter era check (z = 2.983): n_s - 1 = %+.4f" % tilt_numeric(2.983, era_w=0.0), "",
        "PRE-REGISTERED: W1 excluded; W2 with z = 3 exactly excluded if n_s < 1 at > 5 sigma; W2 with fitted eta must show",
        "|alpha_s| < 1e-3 (data must be consistent with zero running); amplitude not predicted."]
txt = "\n".join(out); print(txt); open("iter30_watermark.txt", "w").write(txt + "\n")
