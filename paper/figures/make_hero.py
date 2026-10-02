"""
Hero image of the Grid Universe, in a glowing 3-D style, drawn to match the model:
 - layers = all of space at one moment, each a RANDOM mosaic (not a cubic lattice: a cube grid imprints its axes);
 - a mass dimples every layer it exists in, so its worldline threads straight up through the stack;
 - the fluid between layers (dark matter) pools around the mass;
 - geodesics: an orbiting body's path through space-time is a helix around the mass's worldline (the 'straightest' path in
   the bent grid); light rays crossing a layer bend toward the mass (computed by ray tracing, exaggerated).
Outputs fig0_hero.png (clean) and fig0_hero_labelled.png.
"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import Voronoi
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(11)
GRID, FLUID1, FLUID2, GOLD, LIGHT, WL = "#58aaff", "#3d7bff", "#a35cff", "#ffc94d", "#bfefff", "#e9d9ff"
LAYERS = [0.0, 1.15, 2.3]; A, EPS = 0.42, 0.24; TH = np.radians(28)

def dimple(x, y): return -A/np.sqrt(1 + (x**2 + y**2)/EPS**2)*np.minimum(1, 1.0)
def proj(x, y, z):
    xr = x*np.cos(TH) - y*np.sin(TH); yr = x*np.sin(TH) + y*np.cos(TH)
    s = 1/(1 + 0.16*yr)                                            # mild perspective (far side smaller)
    return xr*s*1.25, z + 0.40*yr*s
def glow(ax, X, Y, color, lw=1.0, layers=((7, .04), (4, .08), (2, .18), (1, .9)), z=2):
    for k, a in layers: ax.plot(X, Y, color=color, lw=lw*k, alpha=a, solid_capstyle="round", zorder=z)
def glow_segments(ax, segs, color, lw=0.8, z=2, base=1.0):
    for k, a in ((6, .035), (3, .08), (1.4, .22), (0.6, .85)):
        ax.add_collection(LineCollection(segs, colors=color, linewidths=lw*k, alpha=a*base, zorder=z, capstyle="round"))

# random mosaic edges in [-1,1]^2 (periodic copies so the border is filled), subdivided so they can bend into the dimple
P = rng.uniform(-1, 1, (230, 2)); Q = np.vstack([P + [dx, dy] for dx in (-2, 0, 2) for dy in (-2, 0, 2)])
vor = Voronoi(Q); edges = []
for rv in vor.ridge_vertices:
    if -1 in rv: continue
    a, b = vor.vertices[rv[0]], vor.vertices[rv[1]]
    if np.all(np.abs(a) <= 1) and np.all(np.abs(b) <= 1): edges.append((a, b))

def layer_segments(z0, amp=1.0):
    segs, nodes = [], []
    for a, b in edges:
        t = np.linspace(0, 1, 9)[:, None]; pts = a + (b - a)*t
        X, Y = proj(pts[:, 0], pts[:, 1], z0 + amp*dimple(pts[:, 0], pts[:, 1]))
        segs += [[(X[i], Y[i]), (X[i + 1], Y[i + 1])] for i in range(len(X) - 1)]
        nodes += [(X[0], Y[0]), (X[-1], Y[-1])]
    # frame
    fr = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]], float)
    return segs, np.array(nodes), proj(fr[:, 0], fr[:, 1], np.full(5, z0))

def render(labels):
    fig = plt.figure(figsize=(12, 9), facecolor="black"); ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor("black")
    ax.set_xlim(-1.95, 1.95); ax.set_ylim(-0.85, 3.05); ax.axis("off")
    # fluid between layers: pools around the mass's worldline (halo-like), colour from blue (low) to violet (high)
    for i in range(len(LAYERS) - 1):
        n = 2600; u = rng.uniform(0, 1, n); r = 0.22*(u/(1 - u + 0.02))**0.75; r = r[r < 1.35]
        ph = rng.uniform(0, 2*np.pi, len(r)); x, y = r*np.cos(ph), r*np.sin(ph); m = (np.abs(x) < 1) & (np.abs(y) < 1)
        x, y = x[m], y[m]; z = LAYERS[i] + rng.uniform(0.12, 0.98, len(x))*(LAYERS[i + 1] - LAYERS[i])
        X, Y = proj(x, y, z); c = FLUID1 if i == 0 else FLUID2
        ax.scatter(X, Y, s=26, color=c, alpha=0.035, lw=0, zorder=1); ax.scatter(X, Y, s=3, color=c, alpha=0.45, lw=0, zorder=1)
    # layers
    for i, z0 in enumerate(LAYERS):
        segs, nodes, (fx, fy) = layer_segments(z0)
        glow_segments(ax, segs, GRID, lw=0.75, z=3 + i)
        ax.scatter(nodes[:, 0], nodes[:, 1], s=4, color="#cfe6ff", alpha=0.55, lw=0, zorder=4 + i)
        glow(ax, fx, fy, GRID, lw=0.9, z=3 + i)
    # mass worldline (straight up through the stack) and the mass in each layer
    zz = np.linspace(LAYERS[0] - 0.35, LAYERS[-1] + 0.45, 200); X, Y = proj(0*zz, 0*zz, zz - A*(np.isin(np.round(zz, 2), LAYERS)))
    glow(ax, *proj(0*zz, 0*zz, zz), WL, lw=1.4, z=9)
    for z0 in LAYERS:
        cx, cy = proj(0, 0, z0 + dimple(0, 0))
        for s, a in ((900, .05), (420, .12), (160, .5)): ax.scatter(cx, cy, s=s, color="#c9b6ff", alpha=a, lw=0, zorder=10)
        ax.scatter(cx, cy, s=110, color="#1a1a22", edgecolor="#e8e4ff", lw=1.2, zorder=11)
        ax.scatter(cx - 0.012, cy + 0.012, s=12, color="white", alpha=0.9, lw=0, zorder=12)
    # geodesic 1: an orbiting body's worldline = helix around the mass's worldline, riding in the dimple
    zt = np.linspace(LAYERS[0] - 0.3, LAYERS[-1] + 0.4, 900); R = 0.42; ph = 2*np.pi*(zt - zt[0])/0.95
    hx, hy = R*np.cos(ph), R*np.sin(ph)
    zh = zt + dimple(R, 0); X, Y = proj(hx, hy, zh); depth = hx*np.sin(TH) + hy*np.cos(TH)   # >0 = far side
    for front, al, zo in ((False, 0.35, 6), (True, 1.0, 13)):
        mk = (depth <= 0) if front else (depth > 0)
        Xm, Ym = np.where(mk, X, np.nan), np.where(mk, Y, np.nan)
        for k, a in ((9, .04), (5, .09), (2.4, .25), (1.2, 1.0)): ax.plot(Xm, Ym, color=GOLD, lw=1.6*k, alpha=a*al, solid_capstyle="round", zorder=zo)
    for z0 in LAYERS:   # where the orbiting body sits in each moment
        j = np.argmin(np.abs(zt - z0)); bx, by = proj(hx[j], hy[j], z0 + dimple(R, 0))
        for s, a in ((260, .08), (90, .3)): ax.scatter(bx, by, s=s, color=GOLD, alpha=a, lw=0, zorder=14)
        ax.scatter(bx, by, s=28, color="#fff3c4", lw=0, zorder=15)
    # geodesic 2: light rays crossing the middle layer, bent toward the mass (ray tracing in n = 1 + k/r, exaggerated)
    def ray(b, k=0.09):
        def f(t, s):
            x, y, vx, vy = s; r = np.hypot(x, y) + 1e-3; n = 1 + k/r; gx, gy = -k*x/r**3, -k*y/r**3
            dot = vx*gx + vy*gy; return [vx/n, vy/n, (gx - dot*vx)/n, (gy - dot*vy)/n]
        s = solve_ivp(f, [0, 2.4], [-1, b, 1, 0], max_step=0.005, rtol=1e-8); m = (np.abs(s.y[0]) <= 1) & (np.abs(s.y[1]) <= 1)
        return s.y[0][m], s.y[1][m]
    zl = LAYERS[1]
    for b in (-0.62, -0.36, 0.36, 0.62):
        x0 = np.linspace(-1, 1, 60); X0, Y0 = proj(x0, 0*x0 + b, zl + dimple(x0, 0*x0 + b)); ax.plot(X0, Y0, ":", color=LIGHT, lw=0.9, alpha=0.45, zorder=8)
        x, y = ray(b, 0.075); X, Y = proj(x, y, zl + dimple(x, y) + 0.004); glow(ax, X, Y, LIGHT, lw=1.3, z=8)
    if labels:
        def lab(xy, txt, col="#dfe9ff", ha="left"):
            ax.text(*xy, txt, color=col, fontsize=11.5, ha=ha, va="center", family="DejaVu Sans", zorder=20, bbox=dict(fc="black", ec="none", alpha=0.75, pad=3))
        lab((1.2, 2.62), "a layer = all of space\nat one moment\n(random mosaic)")
        lab((1.2, 1.78), "fluid between layers\n(dark matter), pooling\naround the mass", "#c8b5ff")
        lab((1.2, 0.9), "light crossing a layer\nbends toward the mass\n(dotted: straight lines)", LIGHT)
        lab((-1.9, 2.75), "the mass's worldline:\nthe same mass, moment\nafter moment", WL)
        lab((-1.9, -0.55), "gold helix: an orbiting body's path\nthrough space-time. It is the\nstraightest path in the bent grid;\nseen one layer at a time, it is a circle", GOLD)
        ax.annotate("", xy=(-1.75, 1.9), xytext=(-1.75, 0.5), arrowprops=dict(arrowstyle="->", color="#9fb3d9", lw=1.4))
        ax.text(-1.8, 1.2, "time", color="#9fb3d9", rotation=90, fontsize=11, va="center", ha="right")
    return fig
for lab, name in ((False, "fig0_hero"), (True, "fig0_hero_labelled")):
    f = render(lab); f.savefig(os.path.join(HERE, name + ".png"), dpi=170, facecolor="black"); plt.close(f); print("wrote", name)
