# Moving the model onto Khronon theory (Blanchet & Skordis 2024, arXiv:2404.06584)

Action: S = c^3/(16 pi G) Int sqrt(-g) [R - 2J(Y) + 2K(Q)] + S_m.  Khronon tau labels a foliation (stack of layers);
n_mu = -(c/Q) grad tau; acceleration A_mu; Y = A.A/c^4; J = Lambda - Y + (2c^2/3a0) Y^{3/2} + ... ; K(Q) = condensate.
Grid picture: foliation = the layers; K(Q) part = fluid between them (dust-like dark matter); J(Y) = MOND response.

Why it avoids AeST's two failures
1. Preferred frame: no separate vector kinetic term (AeST's K_B F^2). In the strong-pull regime J_Y -> 0, so the theory is GR at 1PN:
   alpha1 = -8(alpha-beta)/(1-beta), alpha2 ~ (alpha-beta)(...) vanish with alpha = beta = 0 (paper eqs 6.5-6.6).
   Residual preferred-frame effects scale with the leftover J_Y at Solar-System accelerations (same physics as the Cassini question).
2. Early-universe stability: linear cosmology is generalized dark matter (eq 4.31-4.35): w, c_ad^2 tiny, rest-frame sound speed
   c_s^2 = c_ad^2/(1 + c_ad^2 k^2/(4 pi G a^2 rho (1+w))) >= 0, no anisotropic stress -> no K_B-type tachyon. Only the infrared
   (k < mu ~ 1/22 Mpc) Jeans-like Hamiltonian issue the authors report.
Still open (shared by all MOND theories): Cassini external-field quadrupole vs rotation-curve interpolation.

Implementation: class_aest flag aest_khronon=yes -> GDM equations with k-dependent c_s^2, AeST vector frozen; background and
grid dark energy (beta, exchange) unchanged. Check (khronon_vs_cdm.py): Khronon + Lambda matches LCDM to 0.05% in TT;
insensitive to K_B (no vector). Fits: see fits_summary.txt.
