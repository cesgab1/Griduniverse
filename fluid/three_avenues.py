"""
Three avenues for 'extra mass only in groups/clusters >~1e13.5 Msun, none in galaxies':
Target (fraction of the cosmic share 5.4 x baryons needed at R500; gas_retention_test): 1e12 (MW galaxy, KiDS/SPARC): ~0,
1e13: 0.02, 1e13.5: 0.54, 1e14: 0.71, 1e14.5: 0.57, 1e15: 0.38.
1. Partially locked fluid (drag toward the grid's rest frame, relaxation time tau): a halo keeps fluid only if its pull beats the
   drag on fluid moving with it: g_halo ~ 10 H V_vir > v_rel/tau  ->  trapped if v_rel < k V_esc, k = 10 H tau / sqrt 2 (one number).
   v_rel ~ Maxwellian with halo peculiar speeds. Cost: drag also slows large-scale growth of the fluid when tau H <~ 1.
2. Hot dark component (eV-scale fermion): only deep wells hold it, but Planck+BAO cap its abundance (m_eff < 0.65 eV ->
   Omega_h h^2 < 0.007, i.e. < 6% of the dark matter). A system can hold at most that fraction of the cosmic share.
3. Missing ordinary matter: extra baryons needed = fraction x 5.4 x observed baryons; compare total with the cosmic baryon
   fraction 0.157 of the dynamical mass (a cluster cannot hold more baryons than its share of the universe's).
"""
import numpy as np
from scipy.special import erf
from scipy.optimize import minimize_scalar
exec(open("current_timing.py").read().split("print(f\"{'system'")[0])
target = {12.0: 0.0, 13.0: 0.02, 13.5: 0.54, 14.0: 0.71, 14.5: 0.57, 15.0: 0.38}
fb_obs = {12.0: 0.03, 13.0: 0.097, 13.5: 0.083, 14.0: 0.089, 14.5: 0.116, 15.0: 0.167}
def trapped(Ve, vrms, k):
    x = k*Ve/(vrms/np.sqrt(3)); return erf(x/np.sqrt(2)) - np.sqrt(2/np.pi)*x*np.exp(-x*x/2)
def zform(lM): return 0.9 if lM < 12.5 else (0.7 if lM < 13.8 else 0.5)
def pred1(k): return {lM: trapped(np.sqrt(2)*Vvir(10**lM*1.4, zform(lM)), vpec(zform(lM)), k) for lM in target}
err = lambda lM: 0.15                                                # rough uncertainty on each target value
cost = lambda k: sum(((pred1(k)[lM] - t)/err(lM))**2 for lM, t in target.items())
b = minimize_scalar(cost, bounds=(0.05, 3), method="bounded")
print(f"1. Partially locked fluid: best k = {b.x:.2f} (tau = {b.x*np.sqrt(2)/(10*Hz(0.6))/3.156e16:.2f} Gyr), mismatch chi2 = {b.fun:.1f} (6 targets, 1 parameter)")
for lM, t in target.items(): print(f"      log M {lM}: predicted {pred1(b.x)[lM]:.2f}   needed {t:.2f}")
tau = b.x*np.sqrt(2)/(10*Hz(0.6))
print(f"   cost: drag time tau = {tau/3.156e16:.2f} Gyr vs Hubble time 1/H0 = {1/H0/3.156e16:.1f} Gyr -> tau*H0 = {tau*H0:.2f}; drag at that strength damps the fluid's\n"
      f"   late-time growth on all scales (needs tau*H >> 1 to leave sigma8 alone).")
print("\n2. Hot dark component: maximum fraction of the cosmic share any system can hold = 0.06 (Planck+BAO); needed 0.38-0.71 -> short by 6-12x")
print("\n3. Missing ordinary matter: total baryon fraction of the dynamical mass if the shortfall were baryons")
for lM, t in target.items():
    if lM < 13: continue
    fb_new = fb_obs[lM]*(1 + t*5.4)
    print(f"      log M {lM}: observed {fb_obs[lM]:.3f} -> with missing baryons {fb_new:.2f}  (cosmic {0.157}; {fb_new/0.157:.1f}x the universal ratio)")

# Robustness: X-ray (hydrostatic) masses are thought to be biased low by ~20% (true M500 = M_HSE / 0.8).
print("\nWith a 20% hydrostatic-mass bias (true masses 1.25x higher), the required fraction becomes:")
mond = {13.0: 0.99, 13.5: 0.76, 14.0: 0.66, 14.5: 0.64, 15.0: 0.66}
target_b = {12.0: 0.0}
for lM, m in mond.items():
    mt = m/1.25; fbt = fb_obs[lM]/1.25; target_b[lM] = max(1 - mt, 0)/(5.4*fbt)
print("   " + "  ".join(f"1e{lM}: {t:.2f}" for lM, t in target_b.items()))
cost_b = lambda k: sum(((min(pred1(k)[lM], 1.0) - min(t, 1.0))/0.15)**2 for lM, t in target_b.items())
bb = minimize_scalar(cost_b, bounds=(0.05, 3), method="bounded")
print(f"   partially locked fluid refit: k = {bb.x:.2f}, mismatch chi2 = {bb.fun:.1f}  | predicted: " + "  ".join(f"{pred1(bb.x)[lM]:.2f}" for lM in target_b))
