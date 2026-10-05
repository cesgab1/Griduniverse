import numpy as np
from scipy.integrate import solve_ivp, quad
WR = 2.469e-5*(1+0.2271*3.046)
def model(Om, h, beta):
    Or = WR/h**2; OL = 1 - Om - Or
    if beta == 0: return lambda a: OL
    sol = solve_ivp(lambda x, y: [beta*((Om*np.exp(-3*x)/2 + Or*np.exp(-4*x) - np.exp(y[0]))/(Om*np.exp(-3*x) + Or*np.exp(-4*x) + np.exp(y[0])))],
                    [0, np.log(1e-4)], [np.log(OL)], dense_output=True, rtol=1e-10, atol=1e-12)
    return lambda a: np.exp(sol.sol(np.log(a))[0])
def run(Om, h, beta):
    Or = WR/h**2; rde = model(Om, h, beta)
    E2 = lambda a: Om/a**3 + Or/a**4 + rde(a)
    def rhs(x, y):
        a = np.exp(x); e2 = E2(a); de = 1e-5
        dlnE = (np.log(E2(a*np.exp(de))) - np.log(E2(a*np.exp(-de))))/(4*de)
        Oma = Om/a**3/e2
        return [y[1], -(2 + dlnE)*y[1] + 1.5*Oma*y[0]]
    a0 = 1e-3
    s = solve_ivp(rhs, [np.log(a0), 0], [a0, a0], rtol=1e-9, atol=1e-12)
    D0 = s.y[0, -1]
    age = quad(lambda a: 1/(a*np.sqrt(E2(a))), 1e-8, 1, limit=400)[0]*977.8/(100*h)
    return D0, age
sets = {"Pantheon+": ((0.3016, 0.6845), (0.3115, 0.6733, 0.58)), "DES-Dovekie": ((0.3020, 0.6841), (0.3124, 0.6724, 0.63)),
        "Union3": ((0.3014, 0.6846), (0.3239, 0.6607, 1.19))}
out = []
for k, (L, M) in sets.items():
    DL, aL = run(L[0], L[1], 0); DM, aM = run(*M)
    s8r = DM/DL; S8r = s8r*np.sqrt(M[0]/L[0])
    out.append(f"{k}: H0 constant {100*L[1]:.1f} vs response {100*M[1]:.1f} (local 73.0) | sigma8 ratio {s8r:.3f}, S8 ratio {S8r:.3f} | age constant {aL:.2f} Gyr vs response {aM:.2f} Gyr")
txt = "\n".join(out); print(txt); open("anomalies.txt", "w").write(txt + "\n")
