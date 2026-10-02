"""
Validation of the 3-D Xi solver (pm3d.solve_y): a static MOND configuration must be an equilibrium.
Put a Plummer baryon mass (1e11 Msun) at the box centre at a = 1, give the fluid exactly the MOND phantom density of
that mass (external field 0.025 a0), solve for Xi, and check:
  (1) net force on the fluid, -grad phi + grad Xi, is ~0 compared with gravity (static equilibrium)
  (2) gravity felt by baryons/light, -grad phi, follows MOND: g = nu(gN/a0) gN
Grid effects: 128^3, 47 kpc cells, periodic 6 Mpc box (phantom truncated by the box), so check 100 kpc - 1 Mpc.
NOTE (independent review): the imposed phantom uses a 0.025 a0 external-field estimate but the solver has no external field,
so the ~5% residual net force is a physics mismatch, not a numerical floor; and the 'gravity/MOND' column checks Poisson on the
imposed density only. What this test does validate: the solver converges and the Xi push tracks the fluid's own gravity.
"""
import numpy as np
import pm3d as P
L, N = 6000.0, 128; B = P.Box(L, N); a = 1.0; Mb = 1e11; eps = 5.0
c = (np.arange(N) + 0.5)*L/N; X, Y, Z = np.meshgrid(c, c, c, indexing="ij"); x0 = L/2
r = np.sqrt((X - x0)**2 + (Y - x0)**2 + (Z - x0)**2)
rho_b = 3*Mb/(4*np.pi*eps**3)*(1 + (r/eps)**2)**-2.5
# smooth the baryon mass over a cell so the grid can represent it: Plummer with eps = cell size
eg = B.dx; rho_b = 3*Mb/(4*np.pi*eg**3)*(1 + (r/eg)**2)**-2.5; rho_b *= Mb/(rho_b.sum()*B.dx**3)
rr = np.geomspace(1, 4000, 4000); Mr = Mb*rr**3/(rr**2 + eg**2)**1.5
gN = P.G*Mr/rr**2
# external field 0.025 a0 (1-D estimate) keeps the phantom total (~2e12) well below the box's fluid (7e12), so the zero-mean
# excess never drives the total density negative (an isolated MOND halo in a 6 Mpc box would exceed the cosmic mean)
from scipy.optimize import brentq
gNe = brentq(lambda x: P.nu(x/P.a0)*x - 0.025*P.a0, 1e-12, 10*P.a0)
g = P.nu((gN + gNe)/P.a0)*(gN + gNe) - P.nu(gNe/P.a0)*gNe; Mph = g*rr**2/P.G - Mr
rho_ph = np.gradient(Mph, rr)/(4*np.pi*rr**2)
drho = np.interp(r, rr, rho_ph); drho -= drho.mean()                       # fluid excess = phantom (zero-mean in the box)
drho = P.smooth(B, drho); rho_b = P.smooth(B, rho_b)
y, ds, info = P.solve_y(B, a, drho, np.zeros((N,)*3), tol=1e-4, newton=12, cgmax=80)
phi = B.poisson(4*np.pi*P.G*(drho + rho_b - rho_b.mean())/a)
gphi = B.grad(phi); gy = B.grad(y)
grav = np.sqrt(sum(g_**2 for g_ in gphi)); net = np.sqrt(sum((-gphi[i] - P.cl**2*gy[i])**2 for i in range(3)))
print(f"solver residual {info['res']:.1e}, Newton iterations {info['it']}")
print("  r (kpc) | net force on fluid / gravity | gravity / MOND prediction")
for R in (100, 150, 200, 300, 500, 800, 1200):
    sh = (r > R - B.dx/2) & (r < R + B.dx/2)
    print(f"  {R:6d}  | {np.median(net[sh]/grav[sh]):10.3f}                  | {np.median(grav[sh])/np.interp(R, rr, g):8.3f}")
