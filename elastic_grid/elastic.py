"""Elastic blocks: what tension does space carry, and can gravitational waves be waves of this medium?"""
import numpy as np
c = 2.998e8; rDE = 5.85e-27            # kg/m^3 (measured)
T = rDE*c**2                            # tension = -pressure for w = -1
out = [f"1. Tension of space (dark energy's negative pressure, w = -1): {T:.2e} Pa = {T/101325:.1e} atmospheres."]
out.append("2. A medium whose tension equals its energy density (x c^2) carries transverse waves at EXACTLY c (v^2 = tension/density),")
out.append("   the defining property of relativistic strings/membranes: w = -1 elastic blocks are automatically 'light-speed elastic'.")
# Claim 1 departs from w = -1: d ln rho/d ln a = beta q  ->  w = -1 - beta q/3 ; tension/energy = -w
Om, ODE, beta = 0.31, 0.69, 0.5
for label, q in (("z ~ 1 (decelerating, q ~ +0.2)", 0.2), ("today (q ~ -0.45)", -0.45)):
    w = -1 - beta*q/3; v = np.sqrt(-w)
    out.append(f"3. {label}: w = {w:.3f}, tension/energy = {-w:.3f} -> naive wave speed {v:.4f} c ({100*(v-1):+.1f}%)")
out.append("4. Gravitational waves: speed = c to 1e-15 (GW170817). A 2-4% change is excluded by ~1e13 -> gravitational waves CANNOT be")
out.append("   waves of the dark-energy medium. They must be ripples of the geometry itself (as in GR, iteration 48), and the blocks' own")
out.append("   waves must not propagate: the blocks behave like an INCOMPRESSIBLE, non-ringing elastic medium (no extra 'breathing' wave;")
out.append("   LIGO/Virgo polarisation tests so far consistent with GR's two polarisations only).")
out.append("5. Rate-dependence: rho ~ (stretching rate)^(-1/2): when stretching slows, tension builds; when it speeds up, tension relaxes ->")
out.append("   a self-regulating (thermostat-like) elastic medium; it drives a ~ t^5 in the far future instead of runaway or collapse.")
txt = "\n".join(out); print(txt); open("elastic.txt", "w").write(txt + "\n")
