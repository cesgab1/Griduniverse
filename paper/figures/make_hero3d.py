"""
Cover figure, 3-D version (replaces the rubber-sheet-like dimple drawing).
One moment of space shown in full 3-D. The planes are only a way to draw 3-D space, like the floors of a building; they
are not time. Each plane is a random mosaic. The mass draws the grid in toward itself from EVERY side:
 - the plane through the mass's centre is pinched inward sideways (radially);
 - planes above it sag down toward it, planes below bulge up toward it;
 - the amount of draw-in falls with distance (softened 1/r^2, the shape computed for the random net in Figure 2).
Geodesics: a circular orbit in the central plane; light rays passing above and below bend toward the mass (3-D ray
tracing in the tick-rate field n = 1 + k/r, exaggerated), with dotted straight lines for comparison.
Outputs fig0_grid3d.png (clean) and fig0_grid3d_labelled.png.
"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import Voronoi
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(5)
GRID, GOLD, LIGHT, FLU = "#58aaff", "#ffc94d", "#c4f1ff", "#9b6bff"
PLANES = [-0.84, -0.56, -0.28, 0.0, 0.28, 0.56, 0.84]; C, EPS = 0.05, 0.33
AZ, EL = np.radians(30), np.radians(13)

def pull(p):
    """draw-in toward the mass: displacement r_hat * C/(r^2+eps^2) (monotone: no crossing)."""
    r = np.linalg.norm(p, axis=-1, keepdims=True) + 1e-9; return p - p/r*C/(r**2 + EPS**2)
def proj(p):
    x, y, z = p[..., 0], p[..., 1], p[..., 2]
    xr = x*np.cos(AZ) - y*np.sin(AZ); yr = x*np.sin(AZ) + y*np.cos(AZ)           # yr = depth before tilt
    d = yr*np.cos(EL) - z*np.sin(EL); up = z*np.cos(EL) + yr*np.sin(EL)
    s = 5.0/(5.0 + d); return np.stack([xr*s, up*s], -1), d
def glow_lines(ax, P3, color, width=0.7, zo=3, alpha=1.0):
    P2, d = proj(P3); fade = np.clip(1.05 - 0.35*(d.mean(-1) + 1), 0.35, 1.0)   # far lines dimmer
    for k, a in ((6, .03), (3, .07), (1.5, .2), (0.6, .85)):
        cols = [matplotlib.colors.to_rgba(color, a*alpha*f) for f in fade]
        ax.add_collection(LineCollection(P2, colors=cols, linewidths=width*k, zorder=zo, capstyle="round"))

# random mosaic edges in [-1,1]^2, subdivided so they can bend
P = rng.uniform(-1, 1, (150, 2)); Q = np.vstack([P + [dx, dy] for dx in (-2, 0, 2) for dy in (-2, 0, 2)])
vor = Voronoi(Q); E = []
for rv in vor.ridge_vertices:
    if -1 in rv: continue
    a, b = vor.vertices[rv[0]], vor.vertices[rv[1]]
    if np.all(np.abs(a) <= 1) and np.all(np.abs(b) <= 1): E.append((a, b))
def plane_segments(z0):
    segs = []
    for a, b in E:
        t = np.linspace(0, 1, 8)[:, None]; pts = np.c_[a + (b - a)*t, np.full(8, z0)]
        q = pull(pts); segs += [np.array([q[i], q[i + 1]]) for i in range(7)]
    return np.array(segs)
def ray(p0, v0, k=0.045, T=2.6):
    def f(t, s):
        x = s[:3]; v = s[3:]; r = np.linalg.norm(x) + 1e-3; n = 1 + k/r; g = -k*x/r**3
        return np.r_[v/n, (g - (v @ g)*v)/n]
    s = solve_ivp(f, [0, T], np.r_[p0, v0], max_step=0.004, rtol=1e-9); X = s.y[:3].T
    return X[np.all(np.abs(X[:, :2]) <= 1.0, axis=1)]


# ---- true 3-D random network (no planes, no preferred direction) ----
P3 = np.random.default_rng(9).uniform(-1.15, 1.15, (340, 3)); V3 = Voronoi(P3); NET = []
for rv in V3.ridge_vertices:
    if -1 in rv: continue
    poly = V3.vertices[rv]
    if np.any(np.abs(poly) > 1.0): continue
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]; NET.append((tuple(np.round(a, 9)), tuple(np.round(b, 9))))
NET = np.array(list({tuple(sorted(e)) for e in NET}))          # unique edges
def net_segments():
    segs = []
    for a, b in NET:
        t = np.linspace(0, 1, 6)[:, None]; pts = a + (b - a)*t; q = pull(pts)
        segs += [np.array([q[i], q[i + 1]]) for i in range(5)]
    return np.array(segs)

def render(labels):
    fig = plt.figure(figsize=(12, 9), facecolor="black"); ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor("black")
    ax.set_xlim(-1.62, 1.62); ax.set_ylim(-1.12, 1.12); ax.axis("off")
    # fluid: a round cloud pooled around the mass (dark matter halo), drawn faint
    n = 5000; u = rng.uniform(0, 1, n); r = 0.2*(u/(1 - u + 0.02))**0.8; r = r[(r < 1.1) & (r > 0.05)]
    v = rng.normal(size=(len(r), 3)); v /= np.linalg.norm(v, axis=1, keepdims=True); F = v*r[:, None]
    F = F[np.all(np.abs(F) <= [1, 1, 0.95], axis=1)]; F2, _ = proj(F)
    ax.scatter(*F2.T, s=24, color=FLU, alpha=0.03, lw=0, zorder=1); ax.scatter(*F2.T, s=2.5, color=FLU, alpha=0.35, lw=0, zorder=1)
    # the grid: a random 3-D network (cell edges), every part drawn in toward the mass
    segs = net_segments(); glow_lines(ax, segs, GRID, width=0.55, zo=3, alpha=0.85)
    nodes = np.unique(segs.reshape(-1, 3), axis=0); n2, _ = proj(nodes[::3])
    ax.scatter(*n2.T, s=2.5, color="#d5e8ff", alpha=0.45, lw=0, zorder=4)
    box = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)])
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(box[i] != box[j]) == 1:
                b2, _ = proj(np.array([box[i], box[j]])); ax.plot(*b2.T, color=GRID, lw=0.8, alpha=0.3, zorder=2)
    # geodesic 1: circular orbit in the central plane
    th = np.linspace(0, 2*np.pi, 400); R = 0.5; O = np.c_[R*np.cos(th), R*np.sin(th), 0*th]
    O2, d = proj(O); front = d < 0
    for msk, al, zo in ((~front, 0.35, 5), (front, 1.0, 12)):
        X = np.where(msk, O2[:, 0], np.nan); Y = np.where(msk, O2[:, 1], np.nan)
        for k, a in ((9, .04), (5, .1), (2.4, .28), (1.2, 1)): ax.plot(X, Y, color=GOLD, lw=1.6*k, alpha=a*al, zorder=zo)
    b2, _ = proj(np.array([[R*np.cos(-0.6), R*np.sin(-0.6), 0]]))
    for s, a in ((420, .07), (140, .3), (40, 1)): ax.scatter(*b2.T, s=s, color=GOLD if s > 40 else "#fff3c4", alpha=a, lw=0, zorder=13)
    # geodesic 2: light rays above and below the mass (and one in the plane), bending toward it
    for (y0, z0) in ((0.0, 0.42), (0.0, -0.42)):
        X = ray(np.array([-1.0, y0, z0]), np.array([1.0, 0, 0]), k=0.022); st = np.c_[np.linspace(-1, 1, 60), np.full(60, y0), np.full(60, z0)]
        s2, _ = proj(st); ax.plot(*s2.T, ":", color=LIGHT, lw=1.0, alpha=0.5, zorder=9)
        X2, _ = proj(X)
        for k, a in ((7, .04), (3.5, .1), (1.6, .3), (0.8, 1)): ax.plot(*X2.T, color=LIGHT, lw=1.3*k, alpha=a, zorder=10)
    # the mass
    m2, _ = proj(np.zeros((1, 3)))
    for s, a in ((2600, .04), (1100, .08), (420, .22)): ax.scatter(*m2.T, s=s, color="#c9b6ff", alpha=a, lw=0, zorder=14)
    ax.scatter(*m2.T, s=230, color="#15151d", edgecolor="#efeaff", lw=1.4, zorder=15)
    ax.scatter(m2[0, 0] - 0.018, m2[0, 1] + 0.02, s=22, color="white", alpha=0.9, lw=0, zorder=16)
    if labels:
        def lab(xy, txt, col="#dfe9ff", ha="left"):
            ax.text(*xy, txt, color=col, fontsize=11.5, ha=ha, va="center", zorder=20, bbox=dict(fc="black", ec="none", alpha=0.75, pad=3))
        lab((-1.58, 1.04), "one moment of space, in 3-D:\na random network with no up, down or sideways")
        lab((0.98, 0.8), "above the mass the grid\nis drawn DOWN toward it")
        lab((0.98, 0.02), "level with the mass it is\ndrawn INWARD sideways")
        lab((0.98, -0.8), "below the mass it is\ndrawn UP toward it")
        lab((-1.58, -0.98), "gold: a circular orbit (straightest path\nin the drawn-in grid)\npale blue: light bending toward the mass;\ndotted = straight lines", GOLD)
        lab((-1.58, 0.5), "violet haze: the fluid\n(dark matter) pooled around\nthe mass", "#c8b5ff")
    return fig
for lab, name in ((False, "fig0_grid3d"), (True, "fig0_grid3d_labelled")):
    f = render(lab); f.savefig(os.path.join(HERE, name + ".png"), dpi=170, facecolor="black"); plt.close(f); print("wrote", name)
