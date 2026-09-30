# Can a grid-tension rule make m_e +0.78% at recombination and still pass BBN, quasars, clocks?
import numpy as np
h=0.68; Om=0.31; Or=4.15e-5/h**2; OL=1-Om-Or; beta=0.63
H0yr=h*100/3.0857e19*3.156e7   # H0 in 1/yr
def E(z):   # H/H0 with the beta-law dark energy (rho_DE ∝ adot^-beta), solved iteratively
    a=1/(1+z); e=np.sqrt(Om*a**-3+Or*a**-4+OL)
    for _ in range(50):
        e=np.sqrt(Om*a**-3+Or*a**-4+OL*(a*e)**-beta)
    return e
fr=lambda z: Or*(1+z)/(Or*(1+z)+Om)
zr, zq, zbbn = 1090., 0.886, 4e8
target=0.0078
QSO=4e-7     # 2σ, methanol z=0.886 (Kanekar+2015)
BBN=0.01     # Δm_e/m_e allowed at BBN (Seto+2023)
CLK=1e-16    # |d ln μ/dt| per yr, atomic clocks (order of magnitude)
rules={
 'grid tension = dark energy (β law)': lambda z: np.log(OL*((E(z)/(1+z))**-beta)),
 'tension ∝ expansion rate H (a0 ∝ cH)': lambda z: np.log(E(z)),
 'heat: temperature T': lambda z: np.log(1+z),
}
print(f"{'rule':42s} {'exponent':>9s} {'BBN':>9s} {'quasar z=0.89':>14s} {'clock /yr':>10s}")
for k,f in rules.items():
    p=target/(f(zr)-f(0))
    d=lambda z: p*(f(z)-f(0))
    dz=1e-4; dlnm_dt=p*(f(dz)-f(0))/dz*(-1)*E(0)*H0yr  # dz/dt=-(1+z)H
    print(f"{k:42s} {p:9.2e} {d(zbbn):9.2%} {abs(d(zq)):14.1e} {abs(dlnm_dt):10.1e}")
for n in [0.5,1,2,3]:
    c=target/fr(zr)**n; d=lambda z: c*(fr(z)**n-fr(0)**n)
    dt=c*n*fr(0)**n*E(0)*H0yr*(1-fr(0))  # d fr/dt ≈ -H fr(1-fr)
    print(f"{'radiation strain (f_rad)^'+str(n):42s} {c:9.2e} {c:9.2%} {abs(d(zq)):14.1e} {dt:10.1e}")
print("\nlimits:                                              BBN<1%   quasar<4e-7   clock<1e-16")
print("f_rad at recombination =",round(fr(zr),3), " ratio needed rec->z=0.89:",f"{target/QSO:.0e}")
