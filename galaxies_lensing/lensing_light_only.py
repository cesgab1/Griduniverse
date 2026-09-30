"""
Light bending from light alone: KiDS-1000 weak lensing around isolated galaxies (Brouwer et al. 2021).
x = g_bar: the Newtonian pull of the galaxies' stars (+ cold gas estimate; optionally + hot gas) at each radius, computed by B21 from light.
y = g_obs: the pull measured by how much background galaxies are distorted (lensing), from the ESD: g_obs = 4 G ESD.
Prediction (no fitting to lensing at all): g = g_bar * nu(g_bar/a0), a0 from the blind SPARC rotation-curve test.
Full covariance. Compared with Newton from light alone (no slack, no dark matter).
"""
import numpy as np
G = 4.52e-30; pcm = 3.086e16; D = "/home/claude/kids/"
def load(fn, cov):
    d = np.loadtxt(D + fn); gb, esd, err, bias = d[:, 0], d[:, 1], d[:, 3], d[:, 4]
    gobs = 4*G*esd/bias*pcm
    c = np.loadtxt(D + cov); n = len(gb)
    C = (c[:, 4]/c[:, 6]).reshape(n, n) * (4*G*pcm)**2
    return gb, gobs, C
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
sets = [("isolated galaxies, stars + cold gas (259,000 lenses)", "Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt", "Fig-4-5-C1_RAR-KiDS-isolated_covmatrix.txt"),
        ("same, with estimated hot gas added to g_bar", "Fig-4_RAR-KiDS-isolated_hotgas_Nobins.txt", "Fig-4_RAR-KiDS-isolated_hotgas_covmatrix.txt"),
        ("isolated dwarf galaxies", "Fig-10_RAR-KiDS-isolated-dwarfs_Nobins.txt", "Fig-10_RAR-KiDS-isolated-dwarfs_covmatrix.txt")]
for lab, f, cf in sets:
    gb, go, C = load(f, cf); Ci = np.linalg.inv(C); n = len(gb)
    print(f"\n{lab}: {n} points, g_bar from {gb.min():.1e} to {gb.max():.1e} m/s^2")
    for mlab, pred in (("slack rule, a0 = 0.95e-10 (blind, from SPARC half-split)", gb*nu(gb/0.95e-10)),
                       ("slack rule, a0 = 1.15e-10 (SPARC full fit)", gb*nu(gb/1.15e-10)),
                       ("Newton, light only", gb)):
        r = go - pred; chi = r @ Ci @ r
        ratio = go/pred
        print(f"   {mlab:55s} chi2 = {chi:7.1f} for {n} points | measured/predicted: median {np.median(ratio):.2f}, range {ratio.min():.2f}-{ratio.max():.2f}")
    # where does it deviate: show points
    pred = gb*nu(gb/1.15e-10); e = np.sqrt(np.diag(C))
    print("   g_bar      measured g_obs     predicted     (measured-pred)/error")
    for x, y, p, s in zip(gb, go, pred, e): print(f"   {x:.2e}   {y:.2e} ± {s:.1e}   {p:.2e}   {(y-p)/s:+.1f}")

print("\n=== with the neighbour pull (external field e, Chae et al. 2020 form), one e for the stacked sample ===")
from scipy.optimize import minimize_scalar
def nu_e(y, e):
    A = e*(1+e/2)/(1+e); B = 1+e
    return 0.5 - A/y + np.sqrt((0.5 - A/y)**2 + B/y)
for lab, f, cf in sets[:2]:
    gb, go, C = load(f, cf); Ci = np.linalg.inv(C)
    for a0 in (0.95e-10, 1.15e-10):
        chi = lambda e: (lambda r: r @ Ci @ r)(go - gb*nu_e(gb/a0, e))
        b = minimize_scalar(chi, bounds=(1e-4, 0.5), method="bounded")
        print(f"   {lab[:45]:45s} a0={a0:.2e}: best e = {b.x:.3f}, chi2 {chi(1e-6):.1f} -> {b.fun:.1f} (15 points)")
