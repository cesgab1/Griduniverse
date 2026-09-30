"""
Acceleration-triggered trimming of the ripple (collapse), made testable.
Rule (a CSL-type collapse whose strength depends on how hard the object's centre is accelerated, i.e. its proper acceleration):
   trimming rate per nucleon-mass  lambda(a) = lambda_g * (a/g)^n,   trimming size r_c = 1e-7 m (standard choice)
   For an electron the rate is scaled by (m_e/m_N)^2 (collapse strength grows with mass squared in these models).
Key facts: everything resting on Earth has proper acceleration g (the floor pushes up); free fall (orbit) has ~0.
So on the ground lambda(g) is the ordinary CSL rate: to make measurement pointers collapse it must be >= ~1e-16 /s (GRW value),
and ground experiments cap it (X-ray emission at Gran Sasso: about 5e-12 /s at r_c = 1e-7 m).
Test 1: electrons in storage rings are accelerated ~1e20-1e23 times harder than g for hours. Trimming kicks them randomly;
        the kicks must not spread the beam's energy beyond what is measured (radiation damping keeps it ~0.1%).
Test 2: atoms inside matter (electrons orbiting at ~1e22 m/s^2) -- depends on whether the rule sees internal motion (see note).
"""
import numpy as np
hbar = 1.0546e-34; me = 9.109e-31; mN = 1.6605e-27; c = 2.998e8; g = 9.81; rc = 1e-7; eV = 1.602e-19
rings = {  # name: (energy GeV, bending radius m, damping time s, measured relative energy spread)
    "LEP (CERN, 100 GeV)": (100.0, 3096.0, 0.003, 1.5e-3),
    "SPEAR3 / typical light source (3 GeV)": (3.0, 7.9, 0.004, 1.0e-3),
}
print("Proper acceleration of electrons in rings, and allowed steepness n of the acceleration rule")
for name, (E, R, tau, dE) in rings.items():
    gam = E*1e9*eV/(me*c**2); a = gam**2*c**2/R           # proper acceleration
    p = gam*me*c
    print(f"\n{name}: gamma = {gam:.2e}, proper acceleration = {a:.1e} m/s^2 = {a/g:.1e} g")
    for lam_g, lab in ((1e-16, "weakest that still collapses pointers (1e-16/s)"), (5e-12, "strongest allowed on the ground (5e-12/s)")):
        row = []
        for n in (0, 0.5, 1, 1.5, 2, 3):
            lam = lam_g*(a/g)**n*(me/mN)**2                 # trimming rate for this electron
            D = lam*hbar**2/(2*rc**2)                        # momentum diffusion (kg^2 m^2/s^3), per axis
            sp_ = np.sqrt(D*tau)                              # equilibrium momentum spread against damping
            row.append(f"n={n}: {'OK' if sp_ < dE*p else 'EXCLUDED'} ({sp_/p:.0e})")
        print(f"   lambda(g) = {lab}:  " + ";  ".join(row))
# the window: find max n allowed with lambda_g = 1e-16 for LEP
E, R, tau, dE = rings["LEP (CERN, 100 GeV)"]; gam = E*1e9*eV/(me*c**2); a = gam**2*c**2/R; p = gam*me*c
ns = np.linspace(0, 4, 4001)
ok = [n for n in ns if np.sqrt(1e-16*(a/g)**n*(me/mN)**2*hbar**2/(2*rc**2)*tau) < dE*p]
print(f"\nLEP bound: with the weakest useful ground rate, the rule may grow at most as a^{max(ok):.2f}")
ok2 = [n for n in ns if np.sqrt(5e-12*(a/g)**n*(me/mN)**2*hbar**2/(2*rc**2)*tau) < dE*p]
print(f"           with the strongest allowed ground rate: at most a^{max(ok2):.2f}")
print("\nTest 2 note: if the rule used the electron's internal orbital acceleration (~1e22 m/s^2 in every atom),")
a_atom = (2.19e6)**2/5.29e-11
for n in (0.5, 1, 2):
    lam = 1e-16*(a_atom/g)**n*(me/mN)**2; heat = lam*hbar**2/(2*me*rc**2)   # W per electron (energy gain rate)
    print(f"   n={n}: extra heating {heat*3e26:.1e} W per kg of matter  (Earth's total internal heat ~ 8e-12 W/kg)")
print("\nPrediction unique to the rule: in free fall (orbit) proper acceleration ~0 -> trimming switched OFF.")
print("A space test of a large superposition (e.g. the proposed MAQRO mission) would see NO collapse even if ground tests do.")
