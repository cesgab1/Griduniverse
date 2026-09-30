"""
Can a slack (sparse) grid region explain the CMB Cold Spot?
Facts used (literature): Cold Spot decrement ~ -70 muK (Gaussian-filtered, ~5 deg) up to ~ -150 muK at centre.
Eridanus supervoid (Szapudi et al. 2015): z ~ 0.15-0.2, radius ~ 190-220 Mpc/h, density contrast ~ -0.14 to -0.25.
Standard-model late ISW of that void: about -10 to -30 muK (Szapudi 2015; Kovacs & Garcia-Bellido 2016). So a boost of ~3-7x is needed.
Grid slack rule: pull boosted by nu(g_N/a0); with a background strain e (external field) nu is reduced (Chae et al. 2020 form).
ISW ~ time change of the potential, which the slack rule multiplies by roughly nu (order-of-magnitude treatment).
"""
import numpy as np
a0 = 1.15e-10; Mpc = 3.0857e22; H0 = 67.5e3/Mpc; Om = 0.31; c = 2.998e8
def nu_e(y, e):
    if e == 0: return 0.5 + np.sqrt(0.25 + 1/y)
    A = e*(1 + e/2)/(1 + e); B = 1 + e
    return 0.5 - A/y + np.sqrt((0.5 - A/y)**2 + B/y)
print("Newtonian pull at the edge of cosmic structures (g_N = Om H0^2/2 * |delta| * R):")
for name, R_h, dl in (("Eridanus supervoid", 200, 0.15), ("typical supervoid", 100, 0.3), ("large-scale ISW modes", 500, 0.02)):
    R = R_h/0.675*Mpc; g = Om*H0**2/2*dl*R; y = g/a0
    print(f"   {name:24s} R = {R_h} Mpc/h, |delta| = {dl}: g_N = {g:.1e} m/s^2 = {y:.3f} a0")
    for e, lab in ((0, "no background strain (galaxy rule applied naively)"), (1.0, "background strain = a0"),
                   (9.39e-10/a0, "background strain = cosmic tension g* (8.2 a0)"), (c*H0/a0, "background strain = c*H0 (5.7 a0)")):
        print(f"      boost x{nu_e(y, e):5.2f}   {lab}")
print("\nCold Spot needs ~3-7x the standard ISW of the Eridanus void (-10..-30 muK -> -70..-150 muK).")
print("Whole-sky check: the ISW signal measured by cross-correlating the CMB with galaxy maps has amplitude A ~ 1.0 +- 0.25-0.3")
print("relative to the standard model (Planck), so a boost applying to ALL large-scale structure must be <~1.6x (2 sigma).")
