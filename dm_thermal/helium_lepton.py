"""Nugget formation via a lepton asymmetry: does the helium hint help? (literature check + conversion; NOT pre-registered)"""
from scipy.special import zeta
import numpy as np
xi_e = 0.05                                       # EMPRESS XV (2025): xi_e = 0.05 +0.02 -0.03 (Y_P = 0.2402 +/- 0.0040)
eta_nu = np.pi**2/(12*zeta(3))*(4/11)*(xi_e + xi_e**3/np.pi**2)    # asymmetry per photon
l_e = eta_nu/7.04                                  # per entropy (s = 7.04 n_gamma)
print(f"helium hint xi_e = {xi_e} -> electron-neutrino asymmetry {eta_nu:.4f} per photon = {l_e:.4f} per unit entropy")
print("CMB/BBN limit on TOTAL lepton asymmetry: |l| < 1.2e-2 (Middeldorf-Wygas et al.); flavours can be larger if they cancel")
print("2025 analysis (arXiv:2511.11995): within |l| <~ 1e-2, max baryon chemical potential ~250 MeV; lattice QCD excludes a")
print("  first-order transition below ~450 MeV -> lepton asymmetries CANNOT make the cosmic QCD transition first-order")
print(f"helium-hint asymmetry ({l_e:.4f}) is ~{1.2e-2/l_e:.0f}x below even the maximum allowed, which itself falls short")
