import numpy as np
exec(open("fast_modes.py").read().split("def J(")[0])
def M_of(a, k, KB):
    Q, K, Kp, Kpp, rho, w, H = bg(a); c2a = Kp/(Q*Kpp); Hc = a*H
    def rhs(x):
        d, th, al, E = x
        chi = Q*(a*th/k**2 + al); Sx = KB*E + (2-KB)*chi; Pi = c2a*(d + k**2*Sx/(3*a*a*rho))
        return np.array([-(1+w)*th - 3*Hc*(Pi - w*d), -(1-3*c2a)*Hc*th + k**2*Pi/(1+w), a*E,
                         -Hc*E + a/KB*(Kp*chi - (2-KB)*(Q*Pi/(1+w) + (H+Q)*chi - 3*c2a*H*Q*al))])
    return np.column_stack([rhs(e) for e in np.eye(4)]), Q, a, k, Hc, c2a, rho
for KB in (0.5, 2e-5):
    M, Q, a, k, Hc, c2a, rho = M_of(1e-3, 0.1, KB)
    ev, V = np.linalg.eig(M); i = np.argmax(ev.real); v = V[:, i].real; v = v/np.max(abs(v))
    chi = Q*(a*v[1]/k**2 + v[2])
    print(f"K_B={KB}: growth {ev[i].real/Hc:.2f} aH; eigenvector (delta, theta, alpha, E) = {np.round(v, 6)}; chi = {chi:.3e}; "
          f"pressure source k^2 Sx/(3a^2 rho) relative to delta: {k**2*((2-KB)*chi + KB*v[3])/(3*a*a*rho):.3e}; c2a = {c2a:.1e}")
