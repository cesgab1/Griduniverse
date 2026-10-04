# Pre-registration, iterations 47-49 (written and committed BEFORE any of these runs, Oct 4 2026)

## 47a  The whole history from the grid's own action (no Friedmann equation given to the computer)
Discrete action on the random grid, one term per step between scale factors a_n and a_{n+1} (lapse N_n = time step):
  S = SUM_n N_n abar^3 [ (1/2)<K_ij K^ij - K^2>_grid - rho_m - rho_r - (2/3) C (abar K/3)^(-1/2) ]      (8 pi G = 1)
K_ij read from link-length changes on a random point set (iteration 45 reader). For a homogeneous universe the reader is
exact (every link stretches the same), so this tests the ACTION -> EQUATION step and Claim 1's 'only via K' coupling, not the
reader. The computer finds each N_n by making S stationary (dS/dN_n = 0); the a-equations dS/da_n = 0 are NOT imposed and are
checked afterwards.
 D1 history matches the Friedmann + Claim 1 solution: H(z) within 0.1% (discretisation), crossing z = 0.66 +/- 0.02, w0 = -0.92 +/- 0.01
 D2 the unimposed a-equations hold automatically (residual -> 0 as steps shrink, ~ step^2): the 'second Einstein equation'
    (acceleration) follows from the first plus the action's symmetry. If D2 fails, Claim 1's term is inconsistent.

## 47b  Gravitational waves on a random grid
Ripples h(x, t) on 64,000 random points in a periodic box (spacing 1); the gradient part of the averaged action becomes a
neighbourhood-averaged (Gaussian, width s = 1) operator, normalised analytically (no fitting):
  d^2 h_i/dt^2 = (1/norm) SUM_j exp(-r_ij^2 / 2 s^2) (h_j - h_i)
 G1 wave speed (wavelength 20 spacings) = continuum smoothed prediction c^2 = (1 - e^-x)/x, x = k^2 s^2/2 (c = 0.988) within 1%
 G2 same speed along x, y, z, face and body diagonals: spread < 0.5% (no preferred direction)
 G3 a travelling packet keeps that speed in a real time-stepped run (within 1%)
 G4 consequence for the real grid (cells <= 5.7e-28 m, LIGO wavelengths ~ 3000 km): speed change ~ 1e-69, far below 1e-15.

## 48  Leakage factor c (how much the one frame-dependent field disturbs gravity)
Claim 1 is a function of K (plus the slicing time) -> the VCDM / type-II minimally modified gravity class (BORROWED: two
graviton polarisations only, gravitational-wave speed exactly 1, G_eff = G for sub-horizon modes, static stars = GR/TOV).
 L1 tree level: c = 0 for graviton speed, Solar-System preferred-frame (alpha) and G_eff. Expect this.
 L2 the frame field shifts the K^2 weight for HORIZON-scale perturbations by Delta lambda = Omega_DE/6 (~0.11 today, ~0 at BBN);
    not constrained by the BBN bound; flagged as a horizon-scale effect needing a perturbation calculation.
 L3 loop level: expanding K^(-1/2) around today's K requires the K in Claim 1 to be AVERAGED over large scales (cell-scale K
    fluctuates through zero, where K^(-1/2) is undefined). With averaged K, leakage into gravitons is suppressed to ~0.
    Consequence written now: the graviton-speed bound no longer pins N to the top of its window -> the 'cells right at the
    gamma-ray-burst limit' near-prediction weakens. This is a loss of a test, recorded as such.

## 49  Electron sector (Claim 2)
 E1 Lorentz violation: if the electron mass depends on the tension field only as a scalar (m_e(phi) or m_e(averaged K)),
    there are no direction-dependent electron terms at tree level; LEP bound |c| < 5e-15 satisfied. Expect pass.
 E2 Time variation (the real bite): the switch used in all fits, m_e = 1 + delta/2 [1 + tanh((z - z_s)/30)], leaves a tail.
    Bounds: methanol at z = 0.89, Delta mu/mu = (-1.0 +/- 1.3) e-7 (Bagdonaite et al. 2013 PRL 111, 231101);
    atomic clocks mu_dot/mu ~ 4e-17 /yr. EXPECT: the dz = 30 tail FAILS the z = 0.89 bound for z_s <~ 150 (delta = 0.01), so
    the switch must be sharper (requirement roughly z_s/dz >~ 5). The CMB fit does not care about the width (switch_check).
