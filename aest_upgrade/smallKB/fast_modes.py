"""
Is AeST's cosmology at small K_B physically unstable, or just stiff?
Eigenvalues of the AeST fluid+vector block (delta, theta, alpha, E) with the metric held fixed (fast modes),
same equations as our CLASS patch (Skordis & Zlosnik 2021 eqs 9-12, conformal time), Cosh K(Q), paper K2, Q0, Z0.
Units: Mpc. Background: rho_fld from Q K' - K = 3 rho (CLASS units), H from LCDM-like radiation+matter+DE.
Output: largest real part (growth rate) and largest |imag| (oscillation), in units of the conformal Hubble rate.
"""
import numpy as np
H0 = 68.5/2.998e5; Of, Ob, Or = 0.2512, 0.048, 9.1e-5; OL = 1-Of-Ob-Or
K2, Q0, Z0 = 7.5e3, 0.1, 1e-9
rho0 = Of*H0**2
I0 = 3*rho0/Q0                                 # K' today; K negligible (checked below)
def bg(a):
    s = I0/(2*K2*Z0*a**3); Q = Q0 + Z0*np.arcsinh(s); Kp = I0/a**3; Kpp = 2*K2*np.sqrt(1+s*s)
    K = 2*K2*Z0**2*s*s/(np.sqrt(1+s*s)+1); rho = (Q*Kp-K)/3; w = K/(Q*Kp-K)
    H = H0*np.sqrt((Of+Ob)/a**3 + Or/a**4 + OL); return Q, K, Kp, Kpp, rho, w, H
def J(a, k, KB):
    Q, K, Kp, Kpp, rho, w, H = bg(a); c2a = Kp/(Q*Kpp); Hc = a*H
    M = np.zeros((4, 4))   # rows: d/dtau of (delta, theta, alpha, E); linear in state
    def rhs(x):
        d, th, al, E = x
        chi = Q*(a*th/k**2 + al); Sx = KB*E + (2-KB)*chi; Pi = c2a*(d + k**2*Sx/(3*a*a*rho))
        return np.array([-(1+w)*th - 3*Hc*(Pi - w*d),
                         -(1-3*c2a)*Hc*th + k**2*Pi/(1+w),
                         a*E,
                         -Hc*E + a/KB*(Kp*chi - (2-KB)*(Q*Pi/(1+w) + (H+Q)*chi - 3*c2a*H*Q*al))])
    for j in range(4):
        e = np.zeros(4); e[j] = 1; M[:, j] = rhs(e)
    ev = np.linalg.eigvals(M)
    return ev.real.max()/Hc, np.abs(ev.imag).max()/Hc
print(f"check: K/(Q K') today = {bg(1)[1]/(bg(1)[0]*bg(1)[2]):.1e}")
for KB in (0.5, 0.1, 0.05, 1e-3, 2e-5):
    print(f"\nK_B = {KB}")
    print("   a        k=0.01            k=0.1             k=1   (growth rate / aH , oscillation / aH)")
    for a in (1e-5, 1e-3, 1e-2, 0.1, 1.0):
        row = [J(a, k, KB) for k in (0.01, 0.1, 1.0)]
        print(f"   {a:7.0e}  " + "  ".join(f"{g:+8.2f},{o:9.1e}" for g, o in row))
