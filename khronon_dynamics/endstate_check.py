"""
End-state check for the late Khronon runs (z = 0.7 -> 0, started from the CDM run's z = 0.7 state; before z ~ 0.7 the DBI fluid is
dust-like and Khronon = CDM, see README section 7). For a z = 0 snapshot: solve Xi (with the same external field as the run),
then report around the halo (the centre of mass of the densest region near the box centre), shell by shell:
  enclosed fluid excess, and the radial forces on the fluid: baryon pull g_b, fluid self-gravity g_fl, Xi push,
  net = g_b + g_fl - Xi.  Equilibrium (the MOND state) means net ~ 0; the informative ratio is net / g_b (Xi automatically
  cancels g_fl in the MOND regime, so net/gravity is not informative; independent review, README section 7).
Usage: python3 endstate_check.py <positions.npy> <GEXT in a0>
"""
import numpy as np, sys, os
os.environ["GEXT"] = sys.argv[2] if len(sys.argv) > 2 else "0"
import pm3d as P
L, N = 6000.0, 128; B = P.Box(L, N); m = P.rho_f*L**3/128**3
pos = np.load(sys.argv[1]); drho = P.smooth(B, B.cic(pos, m) - P.rho_f)
c = (np.arange(N) + 0.5)*L/N; X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
sm = P.smooth(B, drho, 2.0); near = (abs(X-L/2) < 800) & (abs(Y-L/2) < 800) & (abs(Z-L/2) < 800)
i = np.unravel_index(np.argmax(np.where(near, sm, -1e30)), sm.shape); x0 = np.array([X[i], Y[i], Z[i]])
dx = [X - x0[0], Y - x0[1], Z - x0[2]]; R = np.sqrt(sum(d**2 for d in dx))
y, ds, info = P.solve_y(B, 1.0, drho, np.zeros((N,)*3), tol=1e-4, newton=20, cgmax=300)
gp = B.grad(B.poisson(4*np.pi*P.G*drho)); gy = B.grad(y)
bx = [X - L/2, Y - L/2, Z - L/2]; rb = np.sqrt(sum(b**2 for b in bx)); fb = P.G*1e11/(rb**2 + 10.0**2)**1.5
rad = lambda F: sum(F[k]*dx[k] for k in range(3))/np.maximum(R, 1e-9)
g_fl = -rad([-g for g in gp]); g_b = -rad([-fb*b for b in bx]); xi = rad([-P.cl**2*g for g in gy])
d = pos - x0; d -= L*np.round(d/L); rp = np.sqrt((d**2).sum(1))
print(f"GEXT = {P.GEXT} a0 | solver res {info['res']:.1e} | halo centre offset {np.round(x0 - L/2)}")
print("  r (kpc) | fluid excess (Msun) | g_b    g_fl   Xi push | net/g_b  (0 = MOND equilibrium; 1 = fluid feels full baryon pull)")
for Rc in (150, 200, 300, 500, 800):
    sh = (R > Rc - B.dx/2) & (R < Rc + B.dx/2)
    exc = m*np.sum(rp < Rc) - P.rho_f*4/3*np.pi*Rc**3
    a_, b_, c_ = np.mean(g_b[sh]), np.mean(g_fl[sh]), np.mean(xi[sh])
    print(f"  {Rc:6d}  | {exc:12.2e}       | {a_:6.2f} {b_:6.2f} {c_:6.2f} | {(a_ + b_ - c_)/a_:6.2f}")
