"""
ITERATION 32 (pre-registration part): Coalesce's inversion -- if dark energy's size is set by the tally of the WHOLE universe
(all of space x all of time), the measured size tells us how big the whole workbook is.
Rule: rho_DE = c_rule / sqrt(N_total), N_total = (comoving volume of ALL space) x integral of a^3 dt over ALL time (Planck units).
Write space as S = V_space / V_visible-today (comoving radius 14.3 Gpc) and let time end at t_end (Branch A: finite book).
Then  S * T(t_end) = (c_rule / rho_DE)^2 / N_visible-past,  with T(t_end) = integral_0^t_end a^3 dt / integral_0^t0 a^3 dt >= 1.
For c_rule = 1 the right side is 4.37. Predictions written BEFORE looking up curvature/topology limits:
 P32a: space is FINITE and small: S < 4.37 (c_rule = 1), i.e. at most ~4.4x the volume we can see.
 P32b: if space is a closed sphere: radius of curvature <= ~14 Gpc -> Omega_k <= -0.1 (expected to be EXCLUDED by Planck+BAO).
 P32c: if space is flat but wraps around (3-torus, side L): L = (S V_vis)^(1/3) between 23 Gpc (S = 1) and 38 Gpc (S = 4.37);
       CMB topology searches would then see repeated patterns if L is close to the last-scattering diameter (~28 Gpc).
 P32d: the larger space is, the sooner time ends (table below). The rule constant c_rule is unknown (O(1)); c_rule = 2 allows
       S up to ~17 (L up to ~60 Gpc), which weakens all of this.
"""
import numpy as np
from scipy.integrate import quad
H0=67.4e3/3.0857e22; Om=0.315; Or=9.1e-5; OL=1-Om-Or; Gyr=3.156e16; c=2.998e8
E=lambda a: np.sqrt(Om/a**3+Or/a**4+OL)
past=quad(lambda a: a**2/(H0*E(a)),1e-10,1)[0]; t0=quad(lambda a: 1/(a*H0*E(a)),1e-10,1)[0]/Gyr
Rv=14.3; Vv=4/3*np.pi*Rv**3
out=["ITERATION 32: whole-workbook inversion (pre-registered before curvature/topology lookup)", "",
     f"visible comoving volume {Vv:.0f} Gpc^3 (radius {Rv} Gpc)", "",
     " c_rule   S (space / visible)   torus side L [Gpc]   closed-sphere radius [Gpc] (Omega_k)   time ends (Gyr from now)"]
for cr in (0.5, 1.0, 2.0):
    rhs=4.37*cr**2
    for S in (1.0, 1.5, 2.0, 3.0, 4.0, rhs*0.999):
        if S>rhs: continue
        T=rhs/S; af=1.0
        while (past+quad(lambda a: a**2/(H0*E(a)),1,af)[0])/past<T: af+=0.002
        tend=quad(lambda a: 1/(a*H0*E(a)),1e-10,af)[0]/Gyr-t0
        L=(S*Vv)**(1/3); Rc=(S*Vv/(2*np.pi**2))**(1/3); Ok=-(c/1e3/67.4/1e3/Rc)**2*1e0
        Ok=-((299792.458/67.4)/1e3/Rc)**2
        out.append(f"  {cr:4.1f}      {S:6.2f}               {L:6.1f}              {Rc:6.1f} ({Ok:+.3f})                    {tend:6.1f}")
    out.append("")
txt="\n".join(out); print(txt); open("iter32_whole_workbook.txt","w").write(txt+"\n")
