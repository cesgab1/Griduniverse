"""
Re-melt check (converged solves only). Take the CDM control's halo at z = 0 (fluid that fell in as dust) and ask whether, under
Khronon with the cosmology-passing DBI law, it is held or pushed out: solve Xi for that density at a = 1 (and, for comparison,
at a = 0.5 and 0.77 with the same comoving density) and measure the shell-averaged radial forces on the fluid around the halo:
gravity (fluid + baryons) and the Xi push. net/gravity ~ 1: the halo is held like CDM; ~ 0: the fluid feels only the baryons
and the excess leaves (re-melt toward the MOND equilibrium).
Also reports the local crossover length 2 pi/mu_eff(local) in the halo, from the solved dust density.
"""
import numpy as np, sys
import pm3d as P
L, N = 6000.0, 128; B = P.Box(L, N); m = P.rho_f*L**3/128**3
pos = np.load("pm3d_cdm_pos.npy"); drho = P.smooth(B, B.cic(pos, m) - P.rho_f)
c = (np.arange(N) + 0.5)*L/N; X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
# halo centre: densest smoothed cell near the box centre
sm = P.smooth(B, drho, 2.0); i = np.unravel_index(np.argmax(np.where((abs(X-L/2)<800)&(abs(Y-L/2)<800)&(abs(Z-L/2)<800), sm, -1e30)), sm.shape)
x0 = np.array([X[i], Y[i], Z[i]]); print("halo centre offset from box centre (kpc):", np.round(x0 - L/2), flush=True)
dx = X - x0[0], Y - x0[1], Z - x0[2]; R = np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)
Mb, eps = 1e11, 10.0
dc = np.array([L/2]*3) - x0
for a in (1.0, 0.77, 0.5):
    y, ds, info = P.solve_y(B, a, drho, np.zeros((N,)*3), tol=1e-3, newton=20, cgmax=300)
    gp = B.grad(B.poisson(4*np.pi*P.G*drho/a)); gy = B.grad(y)
    # radial components (comoving-force units, du/dt); baryons at box centre (physical Plummer force times a)
    bx = [X - L/2, Y - L/2, Z - L/2]; rb = a*np.sqrt(bx[0]**2 + bx[1]**2 + bx[2]**2); fb = a*P.G*Mb*a/(rb**2 + eps**2)**1.5
    grav_r = sum((-gp[k] - fb*bx[k])*dx[k] for k in range(3))/np.maximum(R, 1e-9)
    xi_r = sum((-P.cl**2*gy[k])*dx[k] for k in range(3))/np.maximum(R, 1e-9)
    s1 = P.bg_state(a); _, yp = P.y_of_ds(s1, ds); mu_loc = np.sqrt(P.MU**2/yp)
    print(f"a = {a}: solver res {info['res']:.1e} ({info.get('it')} Newton)")
    print("   r_phys (kpc) | fluid density/mean | gravity (inward) | Xi push (outward) | net/gravity | local crossover (kpc phys)")
    for Rc in (100, 150, 200, 300, 500, 800):
        sh = (R*a > Rc - 0.5*B.dx*a) & (R*a < Rc + 0.5*B.dx*a)
        gr = -np.mean(grav_r[sh]); xp = np.mean(xi_r[sh])
        print(f"   {Rc:8d}     | {np.mean(drho[sh])/P.rho_f + 1:9.1f}         | {gr:10.2f}       | {xp:10.2f}        | {(gr - xp)/gr:7.2f}     | {np.median(2*np.pi/mu_loc[sh]):10.1f}", flush=True)
