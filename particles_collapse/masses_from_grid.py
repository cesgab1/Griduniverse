"""
Can the grid fix particle masses?  Three grid-natural hypotheses, each testable against measured masses.
 M1  'Natural' couplings: every Yukawa coupling is of order 1 at the grid scale. Run down with the two-loop SM equations.
     The top Yukawa is attracted to an infrared quasi-fixed point, so a whole range of grid values gives one top mass.
 M2  Grid-scale seesaw: neutrinos get mass m_nu = y^2 v^2 / M_R with the heavy partner at the grid scale M_R = 1/l_d.
 M3  The hierarchy: what the grid would have to supply to produce the measured Yukawa couplings.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
exec(open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/sm_from_grid.py").read().split("sol = run(inputs(), mu_grid, dense=True)")[0].split('print("\\n=== B & C')[1].split("\n", 1)[1])
v = 246.22
yt_grid = run(inputs(), mu_grid).y[3, -1]
print(f"measured top Yukawa: y_t(m_t) = {inputs()[3]:.4f}  ->  at the grid scale y_t = {yt_grid:.3f}")

print("\n=== M1: natural O(1) Yukawa at the grid scale -> top mass ===")
def run_down(yt_top):
    """shoot: choose y_t(m_t) so that y_t(grid) = yt_top; gauge couplings and lambda from measurement"""
    f = lambda y0: run([*inputs()[:3], y0, inputs()[4]], mu_grid).y[3, -1] - yt_top
    lo, hi = 0.3, 1.25
    for _ in range(60):
        mid = (lo + hi) / 2
        try: val = f(mid)
        except Exception: val = np.inf
        if not np.isfinite(val) or val > 0: hi = mid
        else: lo = mid
    return (lo + hi) / 2
for yg in [0.38, 0.5, 1.0, 2.0, 3.0]:
    y0 = run_down(yg); mt_pred = 173.34 + (y0 - 0.93690) / 0.00556
    print(f"   y_t(grid) = {yg:4.2f}  ->  y_t(m_t) = {y0:.3f}  ->  top pole mass ~ {mt_pred:.0f} GeV")
print("   measured: 172.6 GeV. A 'natural' y_t(grid) = 1-3 gives 200-230 GeV: the fixed point overshoots by ~30-60 GeV.")

print("\n=== M2: grid-scale seesaw ===")
M_R = mu_grid
for name, y in [("top-sized coupling at the grid scale", yt_grid), ("coupling = 1", 1.0)]:
    m_nu = y**2 * v**2 / (2 * M_R) * 1e9          # eV
    print(f"   {name:36s}: heaviest neutrino ~ {m_nu:.1e} eV")
print("   measured (oscillations): the heaviest neutrino is at least sqrt(2.5e-3 eV^2) = 0.05 eV")
print(f"   -> the grid scale is ~{0.05/(1.0**2*v**2/(2*M_R)*1e9):.0e} times too heavy; the seesaw needs M_R ~ {v**2/(2*0.05e-9):.1e} GeV, a scale the grid does not have")

print("\n=== M3: what the hierarchy asks of the grid ===")
masses = {"top": 172.6, "bottom": 4.18, "tau": 1.777, "charm": 1.27, "muon": 0.10566, "strange": 0.0935, "down": 0.00467, "up": 0.00216, "electron": 0.000511}
for n, m in masses.items():
    y = np.sqrt(2) * m / v
    print(f"   {n:9s} Yukawa ~ {y:.1e}  (= 0.22^{np.log(y)/np.log(0.22):.1f})")
print("   The couplings span 5.5 decades and fall roughly on powers of 0.22 (the Cabibbo angle): a Froggatt-Nielsen-type pattern.")
print("   A grid explanation would need an integer 'charge' per particle plus a small grid ratio of 0.22. Nothing in the grid supplies either yet.")
me, mmu, mtau = 0.51099895, 105.6583755, 1776.86
Q = (me + mmu + mtau) / (np.sqrt(me) + np.sqrt(mmu) + np.sqrt(mtau))**2
print(f"\n   Known unexplained regularity (Koide): (m_e+m_mu+m_tau)/(sqrt m_e+sqrt m_mu+sqrt m_tau)^2 = {Q:.6f}  (2/3 = 0.666667)")
