"""Orbits used in the figures (computed, Newtonian limit of the model; units G = 1).
kepler_ellipse : a body around a lone mass -> closed ellipse (Kepler).
rosette        : a body inside the fluid pool (mass + cored halo) -> the orbit does not close; it traces a rosette.
binary_chaos   : a light body near two masses orbiting each other -> irregular, never repeats (chaotic).
"""
import numpy as np
from scipy.integrate import solve_ivp

def _run(acc, x0, v0, T, n=6000):
    f = lambda t, s: np.r_[s[3:], acc(t, s[:3])]
    s = solve_ivp(f, [0, T], np.r_[x0, v0], t_eval=np.linspace(0, T, n), rtol=1e-10, atol=1e-12, method="DOP853")
    return s.y[:3].T, s.t

def point_acc(M=1.0, soft=1e-4):
    return lambda t, x: -M*x/(np.dot(x, x) + soft**2)**1.5

def kepler_ellipse(a=0.55, e=0.6, tilt_deg=32, M=1.0, turns=1.0):
    rp = a*(1 - e); vp = np.sqrt(M*(1 + e)/(a*(1 - e)))
    T = turns*2*np.pi*np.sqrt(a**3/M)
    X, _ = _run(point_acc(M), np.array([rp, 0, 0]), np.array([0, vp, 0]), T)
    c, s = np.cos(np.radians(tilt_deg)), np.sin(np.radians(tilt_deg)); R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return X @ R.T

def halo_acc(M=1.0, Mh=4.0, rc=0.45):
    def acc(t, x):
        r2 = np.dot(x, x) + 1e-8; r = np.sqrt(r2); Menc = M + Mh*r**3/(r2 + rc**2)**1.5
        return -Menc*x/r**3
    return acc

def rosette(r0=0.62, vfac=0.55, tilt_deg=8, T=14.0, **kw):
    acc = halo_acc(**kw); x0 = np.array([r0, 0, 0]); g = np.linalg.norm(acc(0, x0)); vc = np.sqrt(g*r0)
    X, _ = _run(acc, x0, np.array([0, vfac*vc, 0]), T, n=9000)
    c, s = np.cos(np.radians(tilt_deg)), np.sin(np.radians(tilt_deg)); R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return X @ R.T

def binary_chaos(T=60.0, m1=1.0, m2=0.6, sep=0.5, x0=(0.0, 0.62, 0.0), v0=(-0.95, 0.0, 0.08)):
    w = np.sqrt((m1 + m2)/sep**3); r1, r2 = sep*m2/(m1 + m2), sep*m1/(m1 + m2)
    def pos(t):
        return (np.array([r1*np.cos(w*t), r1*np.sin(w*t), 0]), np.array([-r2*np.cos(w*t), -r2*np.sin(w*t), 0]))
    def acc(t, x):
        p1, p2 = pos(t); d1, d2 = x - p1, x - p2
        return -m1*d1/(np.dot(d1, d1) + 1e-4)**1.5 - m2*d2/(np.dot(d2, d2) + 1e-4)**1.5
    X, t = _run(acc, np.array(x0), np.array(v0), T, n=20000)
    P1 = np.array([pos(tt)[0] for tt in t]); P2 = np.array([pos(tt)[1] for tt in t])
    return X, P1, P2
