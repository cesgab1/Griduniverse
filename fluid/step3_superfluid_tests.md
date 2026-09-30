# Step 3: standalone superfluid (Berezhiani-Khoury) against the data

Checked against the published tests before doing our own:

| Test | Result | Source |
|---|---|---|
| Weak lensing out to ~1 Mpc (KiDS-1000, Brouwer+2021) | Fails: reduced chi2 15.0 / 14.9 / 28.7 across three galaxy-mass bins; MOND gets 6.5 with no free parameters | Mistele, McGaugh & Hossenfelder, JCAP 09 (2023) 004, arXiv:2303.08560 |
| Why | The fluid-ripple (phonon) push acts on matter directly; light does not feel it. Rotation (push + Newton) and lensing (Newton only) cannot both look MOND-like | same |
| SPARC rotation curves (169 galaxies) | Fits need star mass-to-light ratios that fall with galaxy size (unnatural); forcing the MOND regime overshoots strong-lensing masses | Mistele & McGaugh, A&A 2022, "Galactic mass-to-light ratios with superfluid dark matter" |
| Solar System | Claimed: superfluidity breaks down near stars. Not tested quantitatively against Cassini in the literature we found | Berezhiani, Famaey & Khoury 2018 |

Verdict: the standalone superfluid is ruled out as-is by lensing, for a structural reason (the push bypasses the geometry light travels through).

## What survives, and where it maps

AeST (already our upgraded base, fitted to CMB+BAO+SN) has the same fluid inside it:
- cosmological scalar sector: 8 pi G rho = Q K'(Q) - K, P = K  -> exactly a P(X) superfluid equation of state (the leftover fluid; behaves as dust early)
- static galaxy sector: deep-MOND term proportional to |Y|^{3/2} -> the same X^{3/2} structure the 1D-channel gas gives (channel_eos.txt)
- the push acts through the metric (vector field = grid direction), so light bends the same way matter moves -> lensing works

So in the grid picture: fluid = AeST's scalar condensate; grid tension/direction = AeST's vector; the fluid pushes matter *by deforming the grid*, not directly.
Open: Cassini (interpolation shape), clusters, value of a0.
