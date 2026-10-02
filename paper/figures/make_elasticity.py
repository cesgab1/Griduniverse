"""
Figure 12: the grid's elasticity, with its limits (min / max), each tied to an equation and to data.
 (a) Stiffness of a tension link: set by Newton's G. Link tension scale F = c^4/G (the Planck force); the stiffness of space
     in Einstein's equations is c^4/(8 pi G). Data fix it today and limit its drift: lunar laser ranging dG/dt / G =
     (7.1 +/- 7.6)e-14 per year (Hofmann & Mueller 2018); big-bang nucleosynthesis |dG/G| < ~10%.
 (b) Cosmic tension (dark energy) over cosmic history, from the Tension-Rate Law rho_DE ~ adot^(-1/2) (BETA = 1/2):
     its maximum is at the Turnover Point; it is smaller in the past and falls to zero in the far future (a ~ t^5).
 (c) Squeeze limit (max compression): the density cap rho_c ~ 0.41 rho_Planck (LQC), compared with real densities.
 (d) Stretch limit: a cell must stay inside its size window (Fig. 10). If cells only stretched with the expansion since
     the Planck era, one Planck cell would now be ~1 mm across, 24 powers of ten over the limit -> the grid must ADD cells as
     space expands (open model requirement).
"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
TXT, MUT, BLUE, GOLD, RED, GRN, VIO = "#e6ecf7", "#93a2bd", "#58aaff", "#ffc94d", "#ff7a7a", "#5fd38d", "#b48cff"
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": TXT, "axes.labelcolor": TXT, "xtick.color": MUT,
                     "ytick.color": MUT, "axes.edgecolor": "#3a4660", "font.size": 9.5, "mathtext.fontset": "dejavusans"})
c, G, hbar = 2.998e8, 6.674e-11, 1.0546e-34
F_P = c**4/G; K_E = c**4/(8*np.pi*G); rho_P = c**5/(hbar*G**2); rho_cap = 0.41*rho_P; lP = np.sqrt(hbar*G/c**3)
H0 = 67.7e3/3.086e22; rho_crit = 3*H0**2/(8*np.pi*G); Om, Or = 0.31, 9.1e-5; rho_DE0 = (1 - Om)*rho_crit
# Tension-Rate Law history (BETA = 1/2): integrate d ln rho_DE / d ln a = BETA q
BETA = 0.5
def rhs(lna, y):
    a = np.exp(lna); r = np.exp(y[0]); rm, rr = Om*a**-3, Or*a**-4; return [BETA*(rm/2 + rr - r)/(rm + rr + r)]
lb = np.linspace(0, -np.log(1101), 3000); lf = np.linspace(0, np.log(30), 1500)
sb = solve_ivp(rhs, [0, lb[-1]], [np.log(1 - Om - Or)], t_eval=lb, rtol=1e-10); sf = solve_ivp(rhs, [0, lf[-1]], [np.log(1 - Om - Or)], t_eval=lf, rtol=1e-10)
L = np.r_[sb.t[::-1], sf.t[1:]]; R = np.exp(np.r_[sb.y[0][::-1], sf.y[0][1:]])/(1 - Om - Or); a = np.exp(L); z = 1/a - 1
imax = np.argmax(R); zmax, Rmax = z[imax], R[imax]; R_rec = R[0]
# stretch: Planck-era cell stretched to today (temperature ratio = stretch factor)
T_P = 1.417e32; stretch = T_P/2.725; l_now = lP*stretch; l_max = 5.7e-28

fig = plt.figure(figsize=(15, 10), facecolor="black")
def panel(rect, title):
    ax = fig.add_axes(rect); ax.set_facecolor("black")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", color=TXT, pad=12); return ax
# (a) stiffness: G's allowed drift, shown as allowed fractional change of link stiffness since the Big Bang
ax = panel([0.05, 0.62, 0.4, 0.27], "(a) Stiffness of a tension link: fixed by G")
ax.set_xscale("log"); ax.set_xlim(1e-4, 1); ax.set_ylim(0, 1); ax.set_yticks([])
llr = (7.1 + 2*7.6)*1e-14*13.8e9
ax.axvspan(1e-4, llr, color="#1d4f2f"); ax.axvspan(llr, 0.1, color="#3a3a1d"); ax.axvspan(0.1, 1, color="#4a1c1c")
ax.text(np.sqrt(1e-4*llr), 0.5, "allowed\ndrift", ha="center", va="center", color="#bff0cf", fontsize=9)
ax.text(np.sqrt(llr*0.1), 0.5, "excluded today\n(lunar laser ranging);\nallowed only in the\nfirst minutes", ha="center", va="center", color="#f2e6a6", fontsize=8.5)
ax.text(0.32, 0.5, "excluded\nalways\n(nucleo-\nsynthesis)", ha="center", va="center", color="#ffc9c9", fontsize=8.5)
ax.set_xlabel("allowed change of link stiffness over 13.8 billion years  |ΔG/G|")
fig.text(0.05, 0.505, rf"link tension scale  $F = c^4/G = {F_P:.2e}$ N  (Planck force);   stiffness of space  $c^4/8\pi G = {K_E:.1e}$ N"
        "\n" r"drift: $\dot G/G = (7.1 \pm 7.6)\times10^{-14}\,\mathrm{yr}^{-1}$  $\Rightarrow$  " rf"$|\Delta G/G| < {llr:.1e}$ (2σ) over the age of the universe",
        fontsize=9, color=TXT, va="top")
# (b) cosmic tension over history
ax = panel([0.55, 0.62, 0.42, 0.27], "(b) Cosmic tension over history (Tension-Rate Law)")
ax.plot(1 + z, R, color=BLUE, lw=2.2, label=r"grid: $\rho_{DE}\propto \dot a^{-1/2}$"); ax.axhline(1, color=MUT, ls="--", lw=1.3, label="Einstein's constant Λ")
ax.set_xscale("log"); ax.invert_xaxis(); ax.set_xlabel("1 + redshift   (past ← → future)"); ax.set_ylabel("tension / today's")
ax.plot(1 + zmax, Rmax, "o", color=GOLD, ms=8, mec="black"); ax.annotate(f"MAX: {Rmax:.3f}× today's\nat the Turnover Point (z = {zmax:.2f})", (1 + zmax, Rmax), (300, 1.02), color=GOLD, fontsize=9,
            arrowprops=dict(arrowstyle="-", color=GOLD))
ax.annotate(f"at the CMB (z = 1100):\n{R_rec:.2f}× today's", (1101, R_rec), (300, 0.62), color=TXT, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUT))
ax.text(0.05, 0.62, "MIN: falls toward zero in\nthe far future (a ~ t⁵)", color=TXT, fontsize=8.5, transform=ax.transAxes)
ax.legend(frameon=False, fontsize=8.5, loc="lower center", labelcolor=TXT)
fig.text(0.55, 0.505, rf"today: $\rho_{{DE}}c^2 = {rho_DE0*c**2:.1e}$ J/m³ (= tension per area per length, Pa);  "
        "\n" r"law: $d\ln\rho_{DE}/d\ln a = \frac{1}{2}\,q$,  $q = -\ddot a\,a/\dot a^2$;  MAX where $q = 0$ (the Turnover Point)", fontsize=9, color=TXT, va="top")
# (c) squeeze limit
ax = panel([0.05, 0.16, 0.4, 0.24], "(c) Squeeze limit: the density cap (MAX compression)")
ax.set_xscale("log"); ax.set_xlim(1e-30, 1e100); ax.set_ylim(0, 1); ax.set_yticks([])
for x, lab, col, y in ((rho_crit, "average\nuniverse", MUT, 0.75), (1e3, "water", MUT, 0.75), (2.3e17, "atomic\nnucleus", MUT, 0.75),
                       (1e18, "neutron-star\ncore", MUT, 0.3), (rho_cap, f"DENSITY CAP\n{rho_cap:.1e} kg/m³", RED, 0.75)):
    ax.axvline(x, color=col, lw=1.6 if col == RED else 0.9); ax.text(x*3, y, lab, color=col if col == RED else TXT, fontsize=8.5, va="center")
ax.axvspan(rho_cap, 1e100, color="#4a1c1c"); ax.set_xlabel("density (kg/m³)")
fig.text(0.05, 0.075, r"$\rho_{cap} \approx 0.41\,\rho_P$,  $\rho_P = c^5/\hbar G^2$;  collapse obeys $H^2 = \frac{8\pi G}{3}\rho\,(1 - \rho/\rho_{cap})$:"
        "\nat the cap the squeeze turns into a rebound (Bounce, Planck star).\nAssumed from loop quantum cosmology, not measured.",
        fontsize=9, color=TXT, va="top")
# (d) stretch limit
ax = panel([0.55, 0.16, 0.42, 0.24], "(d) Stretch limit: cells cannot simply stretch")
ax.set_xscale("log"); ax.set_xlim(1e-37, 1e0); ax.set_ylim(0, 1); ax.set_yticks([])
ax.axvspan(lP, l_max, color="#1d4f2f"); ax.axvspan(l_max, 1, color="#4a1c1c")
ax.text(np.sqrt(lP*l_max), 0.75, "allowed\ncell size", ha="center", color="#bff0cf", fontsize=9)
ax.axvline(l_now, color=GOLD, lw=1.8); ax.text(l_now/3e8, 0.25, f"one Planck cell stretched\nwith the universe since the\nPlanck era: {l_now:.0e} m", color=GOLD, fontsize=8.5)
ax.set_xlabel("size of one cell today (m)")
fig.text(0.55, 0.075, rf"stretch since the Planck era $= T_P/T_{{CMB}} = {stretch:.1e}$;  allowed stretch per cell $\leq \ell_{{max}}/\ell_P = {l_max/lP:.1e}$"
        "\n" r"$\Rightarrow$ the grid must ADD cells as space expands (OPEN requirement of the model).", fontsize=9, color=TXT, va="top")
fig.text(0.01, 0.985, "Figure 12. The grid's elasticity and its limits: how stiff, how tense, how far it can be squeezed and stretched",
         fontsize=13.5, fontweight="bold", color=TXT, va="top")
for ext in ("png", "pdf"): fig.savefig(os.path.join(HERE, f"fig12_elasticity.{ext}"), dpi=160, facecolor="black")
print(f"F_P {F_P:.3e} N, stiffness {K_E:.2e} N, max tension {Rmax:.3f} at z {zmax:.2f}, at CMB {R_rec:.3f}, rho_cap {rho_cap:.2e}, "
      f"stretch {stretch:.2e}, Planck cell today {l_now:.1e} m, LLR |dG/G| {llr:.1e}")
