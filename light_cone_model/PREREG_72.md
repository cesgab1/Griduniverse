# Pre-registration, iteration 72: reverse calculation from today's measured dark energy (Oct 4 2026, written BEFORE running)

Rename: capital C (dark-energy strength) is now D0 ('dark-energy dial'), to avoid confusion with c (speed of light).

Inputs (all measured, none tuned): Omega_DE = 0.685, Omega_k = 0.0023, Omega_r = 9.1e-5, Omega_m = rest, H0 = 67.4 km/s/Mpc,
fading law beta = 1/2 in the one equation  E^2 = Om(1+z)^3 + Or(1+z)^4 + Ok(1+z)^2 + ODE (E/(1+z))^(-1/2).
Run it backward (and forward for the event horizon). Compute: crossing (rho_DE = rho_m), start of acceleration, ages.

Clue checks and expected outcomes:
R1 (idea 2, 'matter hits the quantum ceiling'): rho_m x horizon^2 at the crossing vs the far-future value 0.1865 c^4/G.
    Two horizon choices (event horizon, Hubble radius) = 2 free choices. PASS if within 30% for the event horizon.
    Also compute the SAME number for a plain constant (Lambda) -- if Lambda gives the same, the hit carries no information.
    Expectation: FAIL (no mechanism forces it); if it passes, expect Lambda to pass too.
R2 (idea 3, 'black holes start the clock'): time of crossing minus time of the star-formation peak (z = 1.86, Madau & Dickinson 2014).
    Clue only (no predicted delay exists). Expectation: delay 2.5-3.5 Gyr.
R3 D0 in physical terms: D0 in units of the Planck density, and its energy scale D0^(1/4). Expectation ~3.3 meV.
    WARNING (look-elsewhere): any closeness to a particle mass (e.g. our derived lightest neutrino 2.24 meV) is flagged, NOT counted,
    because there are dozens of masses to compare with.
R4 Does the empty light-cone state (aH = c), where D0 is defined, ever occur in the past or future? Expectation: never (D0 is a
    reference value, not a moment in history).
R5 (idea 4, information count): horizon entropy at the crossing vs (a) CMB photon+neutrino entropy inside the horizon,
    (b) black-hole entropy today ~3e104 (Egan & Lineweaver 2010). PASS only if a ratio is within x3 of 1 with no powers taken.
    Expectation: FAIL (horizon ~1e122 is far larger than both).
