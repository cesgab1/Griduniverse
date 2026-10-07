# PREREG 124: does the lock leak? (exclusion as a write-lock on the bit plane)
Committed before collecting data and before code.

Idea (Coalesce): particles read/write the bit plane; fermions take an exclusive write lock, bosons share.
If the lock is a mechanism (not an exact law), it may leak. Leak should depend on how a particle's reading footprint
compares with the grid cell.

Model L: probability that one interaction ignores the lock  P = (l / lambda)^k
  l = grid cell size (metres, unknown); lambda = hbar/(m c) = reading footprint of capacity m (measured masses only);
  k = 1 or 2 (free choice; both reported).
Data (published, to collect): Pauli-violation limits for electrons (VIP-2 type: new electrons into filled atoms),
for nucleons (Borexino / Super-K / other nuclear transitions), and for photons (Bose-symmetry violation: a boson
being forced to 'lock').
Output: upper limit on l for each particle and k.
Expectations (written now): electrons, k = 1 -> l far below Planck size (lock cannot leak linearly);
k = 2 -> electron data weak; nucleon data (much smaller limits) exclude even k = 2 at ~Planck-size cells.
Caveat recorded in advance: standard theory forbids a fixed particle changing its family (Messiah-Greenberg), so
nucleon limits test a particular kind of leak; reported, with the electron result as headline.
Free choices: k (2 values), which limit per particle (strongest published, plus one alternative).
