"""Figure 11: three kinds of orbit in the pulled grid (computed, orbits.py), seen from above."""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from orbits import kepler_ellipse, rosette, binary_chaos, halo_acc
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(2)
GOLD, ROS, RED, VIO, TXT, MUT = "#ffc94d", "#4dffc8", "#ff7a7a", "#9b6bff", "#e6ecf7", "#93a2bd"
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": TXT})
def glow(ax, X, Y, col, lw=1.2, al=1.0):
    for k, a in ((8, .035), (4, .09), (2, .25), (1, 1)): ax.plot(X, Y, color=col, lw=lw*k, alpha=a*al, solid_capstyle="round")
def mass(ax, x, y, s=160):
    for ss, a in ((s*6, .05), (s*2.5, .15)): ax.scatter(x, y, s=ss, color="#c9b6ff", alpha=a, lw=0)
    ax.scatter(x, y, s=s, color="#15151d", edgecolor="#efeaff", lw=1.3, zorder=5)
fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), facecolor="black")
for ax in axs: ax.set_facecolor("black"); ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.85, 1.3)
# (a) ellipse around a lone mass: 5 turns overlap exactly
E = kepler_ellipse(tilt_deg=0, turns=5); ax = axs[0]; E = E*1.4; glow(ax, E[:, 0] + 0.46, E[:, 1], GOLD, 1.3); mass(ax, 0.46, 0)
ax.set_title("(a) Elliptical orbit around a lone mass", color=TXT, fontsize=12.5, fontweight="bold", loc="left")
ax.text(-1.25, -1.8, "five laps drawn: they land exactly on top of each other.\nThe mass sits at one focus; the body speeds up when\nclose and slows when far (Kepler's laws).", fontsize=9.5, color=TXT, va="bottom")
# (b) rosette inside a fluid pool
ax = axs[1]; n = 5000; u = rng.uniform(0, 1, n); r = 0.45*(u/(1 - u + 0.02))**0.6; r = r[r < 1.2]; th = rng.uniform(0, 2*np.pi, len(r))
ax.scatter(r*np.cos(th), r*np.sin(th), s=2.2, color=VIO, alpha=0.35, lw=0); ax.scatter(r*np.cos(th), r*np.sin(th), s=20, color=VIO, alpha=0.025, lw=0)
R = rosette(tilt_deg=0, T=14.0)*1.5; glow(ax, R[:, 0], R[:, 1], ROS, 0.9); mass(ax, 0, 0, 110)
ax.set_title("(b) Irregular orbit inside the fluid pool", color=TXT, fontsize=12.5, fontweight="bold", loc="left")
ax.text(-1.25, -1.8, "the pull grows as the body dips deeper into the fluid,\nso each lap swings round further: the orbit never\ncloses and fills a ring (how stars move in galaxies).", fontsize=9.5, color=TXT, va="bottom")
# (c) chaotic orbit near two masses orbiting each other
X, P1, P2 = binary_chaos(T=60, x0=(0, 1.05, 0), v0=(1.0, 0, 0.05)); ax = axs[2]
glow(ax, X[:, 0], X[:, 1], RED, 0.7, 0.9); ax.plot(P1[:, 0], P1[:, 1], color=MUT, lw=0.8, alpha=0.6); ax.plot(P2[:, 0], P2[:, 1], color=MUT, lw=0.8, alpha=0.6)
mass(ax, *P1[-1, :2], 130); mass(ax, *P2[-1, :2], 90)
Y, _, _ = binary_chaos(T=60, x0=(0, 1.05 + 1e-6, 0), v0=(1.0, 0, 0.05)); grow = np.linalg.norm(X - Y, axis=1).max()/1e-6
ax.set_title("(c) Irregular orbit near two masses", color=TXT, fontsize=12.5, fontweight="bold", loc="left")
ax.text(-1.25, -1.8, f"two masses circle each other (grey); a light body is\ntugged by both and never repeats. Shift its start by\none part in a million and the gap grows ~{grow:,.0f}-fold.", fontsize=9.5, color=TXT, va="bottom")
fig.text(0.005, 0.985, "Figure 11. Three kinds of orbit in the pulled grid (computed; seen from above)", fontsize=13.5, fontweight="bold", color=TXT, va="top")
fig.text(0.01, -0.02, r"Equations:  $\ddot{\mathbf{x}} = -\nabla\phi$.   (a) $\phi = -GM/r$: $r = a(1-e^2)/(1+e\cos\theta)$, $T^2 = 4\pi^2 a^3/GM$.   "
         r"(b) $M(r) = M + M_h r^3/(r^2+r_c^2)^{3/2}$, $\ddot{\mathbf{x}} = -G M(r)\,\hat r/r^2$.   (c) $\phi = -\sum_k G m_k/|\mathbf{x}-\mathbf{x}_k(t)|$; gap $\delta(t)\approx\delta_0 e^{\lambda t}$.", fontsize=10.5, color=TXT)
for ext in ("png", "pdf"): fig.savefig(os.path.join(HERE, f"fig11_orbits.{ext}"), dpi=160, facecolor="black", bbox_inches="tight")
print("growth", grow)
