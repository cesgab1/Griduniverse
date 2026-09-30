# a0 from riptide exchange counting, vs local a0 and MUSE-DARK III (2026) a0(z) = 1.0e-10 + 1.59e-10 z
import numpy as np
c=2.998e8; Mpc=3.0857e22; H0=67.7e3/Mpc; Om=0.311; OL=1-Om; beta=0.5; G=6.674e-11
def E(z):
    a=1/(1+z); e=np.sqrt(Om*a**-3+OL)
    for _ in range(200): e=np.sqrt(Om*a**-3+OL*(a*e)**-beta)
    return e
rDE=lambda z: (E(z)/(1+z))**-beta          # rho_DE(z)/rho_DE0
cands={
 'tension energy: g²/8πG = ρ_DE c²  (earlier)':      (lambda z: c*np.sqrt(8*np.pi*G*OL*3*H0**2/(8*np.pi*G)*rDE(z))),
 'exchange bias, straight crossings: pN = 1 → 2cH': (lambda z: 2*c*H0*E(z)),
 'exchange bias, loop circulation: → cH':            (lambda z: c*H0*E(z)),
}
local=1.15e-10; muse=lambda z: 1.0e-10+1.59e-10*z
print(f"local a0 (SPARC fit) {local:.2e};  MUSE-DARK z=1: {muse(1):.2e} (ratio {muse(1)/1.0e-10:.2f} ± ~0.15)\n")
for k,f in cands.items():
    print(f"{k:48s} a0(0) = {f(0):.2e} ({f(0)/local:4.1f}× local)   a0(1)/a0(0) = {f(1)/f(0):.2f}   a0(1.4)/a0(0) = {f(1.4)/f(0):.2f}")
print(f"\nMUSE-DARK ratios: z=1 {muse(1)/muse(0):.2f}, z=1.4 {muse(1.4)/muse(0):.2f}")
