"""CAMB with dark energy following the throat-flow law  d ln rho_DE / d ln a = beta * q  (w = -1 - beta q / 3)."""
import numpy as np
from scipy.integrate import solve_ivp
from cobaya.theories.camb import CAMB

class CAMBLaw(CAMB):
    beta: float = 0.5
    def set(self, params_values_dict, state):
        pars = super().set(params_values_dict, state)
        if pars is None:
            return pars
        h = pars.H0 / 100.0
        om = (pars.ombh2 + pars.omch2 + pars.omnuh2) / h**2
        orad = 4.18e-5 / h**2
        ol = 1 - om - orad
        lna_grid = np.linspace(0.0, np.log(1e-5), 600)
        def rhs(lna, y):
            a = np.exp(lna); rde = np.exp(y[0]); rm, rr = om * a**-3, orad * a**-4
            return [self.beta * (rm / 2 + rr - rde) / (rm + rr + rde)]
        sol = solve_ivp(rhs, [0.0, lna_grid[-1]], [np.log(ol)], t_eval=lna_grid, rtol=1e-8, atol=1e-10)
        a = np.exp(sol.t); rde = np.exp(sol.y[0]); rm, rr = om * a**-3, orad * a**-4
        q = (rm / 2 + rr - rde) / (rm + rr + rde)
        w = -1 - self.beta * q / 3
        pars.set_dark_energy_w_a(a[::-1], w[::-1], dark_energy_model="ppf")
        return pars
