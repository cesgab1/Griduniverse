"""
Neutrino oscillation in the River-Compactness framework, and what the data demand.
 A. Oscillation as grid-step counting: each mass state accumulates phase m c^2 tau / hbar, tau = proper time = grid steps along its path.
    For ultra-relativistic neutrinos tau_i = L m_i c / E, so the phase difference is dm^2 c^3 L / (2 hbar E): the standard formula.
    Checked against the measured oscillation experiments (reactor, long-baseline), and the gravity (compactness) correction on Earth.
 B. Reverse-engineer the neutrino mass matrix in the flavour basis from the measured angles and splittings: what texture is required?
 C. Test the missing ingredient: heavy right-handed neutrinos (grid-allowed: SM singlets add no anomaly) with Z_N charges -> seesaw.
Global-fit oscillation parameters (NuFIT-like, normal ordering): sin^2 th12 0.307, sin^2 th13 0.0220, sin^2 th23 0.50 (octant open),
delta_CP ~ 1.2 pi, dm21 = 7.42e-5 eV^2, dm31 = 2.51e-3 eV^2.
"""
import numpy as np, itertools, json
rng = np.random.default_rng(11)
s12, s13, s23, dcp, dm21, dm31 = 0.307, 0.0220, 0.50, 1.2*np.pi, 7.42e-5, 2.51e-3
def pmns(s12, s13, s23, d, a21=0.0, a31=0.0):
    c12, c13, c23 = np.sqrt(1-s12), np.sqrt(1-s13), np.sqrt(1-s23); S12, S13, S23 = np.sqrt(s12), np.sqrt(s13), np.sqrt(s23)
    U = np.array([[c12*c13, S12*c13, S13*np.exp(-1j*d)],
                  [-S12*c23 - c12*S23*S13*np.exp(1j*d), c12*c23 - S12*S23*S13*np.exp(1j*d), S23*c13],
                  [S12*S23 - c12*c23*S13*np.exp(1j*d), -c12*S23 - S12*c23*S13*np.exp(1j*d), c23*c13]])
    return U @ np.diag([1, np.exp(1j*a21/2), np.exp(1j*a31/2)])
U = pmns(s12, s13, s23, dcp)

print("=== A. Oscillation = difference in grid-step counts between mass states ===")
hbarc = 1.973269804e-7   # eV m
def prob(a, b, L_km, E_GeV, anti=False):
    m2 = np.array([0, dm21, dm31]); Uc = U.conj() if anti else U
    ph = np.exp(-1j * m2 * (L_km*1e3) / (2 * E_GeV*1e9 * hbarc))     # phase_i = m_i^2 c^3 L / (2 hbar E)
    return abs(np.sum(Uc[b] * ph * Uc[a].conj()))**2
for name, a, b, L, E, anti, meas in [("Daya Bay far detector, reactor anti-nu_e survival", 0, 0, 1.65, 0.004, True, "~0.92 (deficit ~8%)"),
                                     ("KamLAND, ~180 km, 4 MeV anti-nu_e survival", 0, 0, 180, 0.004, True, "~0.6 (averaged)"),
                                     ("T2K nu_mu survival at the oscillation dip (295 km, 0.6 GeV)", 1, 1, 295, 0.6, False, "~0.0-0.05"),
                                     ("T2K nu_mu -> nu_e appearance (295 km, 0.6 GeV)", 1, 0, 295, 0.6, False, "~0.05")]:
    if "KamLAND" in name:
        Pav = np.mean([prob(a, b, LL, EE, anti) for LL in np.linspace(140, 215, 40) for EE in np.linspace(0.0025, 0.007, 40)])
        print(f"   {name}: P = {Pav:.3f} averaged over 140-215 km and 2.5-7 MeV   (observed {meas})"); continue
    print(f"   {name}: P = {prob(a, b, L, E, anti):.3f}   (observed {meas})")
C_earth = 6.6743e-11 * 5.972e24 / (6.371e6 * 2.99792458e8**2)
print(f"   Earth's gravity (compactness C = {C_earth:.1e}) shifts every phase by the same factor 1 + C: effect on oscillations ~ {C_earth:.0e} -> invisible")
print("   Grid discreteness: the dispersion correction is the same for all three mass states (flavour-blind), so it cancels in the phase differences")

print("\n=== B. What mass matrix do the data require? |M_ab| in meV, flavour basis (e, mu, tau) ===")
for m1 in (0.0, 0.01, 0.03):
    m = np.array([m1, np.sqrt(m1**2 + dm21), np.sqrt(m1**2 + dm31)])
    for a21, a31 in [(0, 0), (np.pi, 0)]:
        Ua = pmns(s12, s13, s23, dcp, a21, a31)
        M = Ua.conj() @ np.diag(m) @ Ua.conj().T
        A = np.abs(M) * 1e3
        eps = 0.2; nexp = np.log(np.abs(M) / np.abs(M).max()) / np.log(eps)
        print(f"   m1 = {m1*1e3:4.0f} meV, Majorana phase {a21/np.pi:.0f}pi: |M| =", np.round(A, 1).tolist(), " ~ eps^", np.round(nexp, 1).tolist())
print("   -> the mu-tau block is uniformly large (entries equal within ~30%) while the e-row is ~eps smaller, AND the mu-tau block is")
print("      nearly rank 1 (its determinant is suppressed: m2/m3 = 0.17). Equal lepton-doublet charges for mu and tau give the first feature;")
print("      the second needs a rank-reducing mechanism, which a Weinberg operator with O(1) coefficients does not supply.")
