"""Energy-conserving version of the early stronger-G phase: an extra component X with
rho_X(a) = (G(a)/G0 - 1) * rho_rest(a), G(a)/G0 = 1 + delta/(1+(a/a_t)^4), carried together with dark energy
(Lambda or the throat-flow law) as one component. Background expansion is exactly that of the boosted-G universe;
X's perturbations behave like a fluid with sound speed^2 = CS2 (1 = stiff scalar field, 1/3 = radiation-like; PPF for the law).
Helium from BBN with the boost (equivalent dN_eff = 6.14 delta)."""
import os, ctypes, numpy as np
from scipy.integrate import solve_ivp
import camb.bbn as bbn
from camb.baseconfig import camblib
from cobaya.theories.camb import CAMB
camblib.set_gboost.argtypes = [ctypes.c_double, ctypes.c_double]; camblib.set_gboost(0.0, 1e-4)
DELTA = float(os.environ.get("GDELTA", 0)); AT = float(os.environ.get("GAT", 1e-4)); CS2 = float(os.environ.get("GCS2", 1.0))
_bbn = bbn.get_predictor()
def table(pars, beta):
    h = pars.H0/100; om = Om = (pars.ombh2 + pars.omch2 + pars.omnuh2)/h**2; orad = 4.18e-5/h**2; ol = 1 - om - orad
    lna = np.linspace(0, np.log(1e-9), 4000); a = np.exp(lna); rrest = om*a**-3 + orad*a**-4
    if beta is None: rde = np.full_like(a, ol)
    else:
        f = lambda x, y: [beta*(om*np.exp(-3*x)/2 + orad*np.exp(-4*x) - np.exp(y[0]))/(om*np.exp(-3*x) + orad*np.exp(-4*x) + np.exp(y[0]))]
        rde = np.exp(solve_ivp(f, [0, lna[-1]], [np.log(ol)], t_eval=lna, rtol=1e-9, atol=1e-11).y[0])
    if os.environ.get("GFORM") == "recoil":   # stored stretch energy -> released as fast recoil motion (rho ~ a^-6) after a_t
        rrest_t = Om*AT**-3 + orad*AT**-4
        rt = rde + DELTA/(1/rrest + (a/AT)**6/rrest_t)
    elif os.environ.get("GFORM") == "waves":  # released as grid waves (radiation, rho ~ a^-4) after a_t
        rrest_t = Om*AT**-3 + orad*AT**-4
        rt = rde + DELTA/(1/rrest + (a/AT)**4/rrest_t)
    else:
        rt = rde + DELTA/(1 + (a/AT)**4)*rrest
    w = -1 - np.gradient(np.log(rt), lna)/3
    if beta is None: w = np.maximum(w, -1 + 1e-9)
    return a[::-1], w[::-1]
class _Base(CAMB):
    beta = None
    def set(self, params_values_dict, state):
        pars = super().set(params_values_dict, state)
        if pars is None: return pars
        pars.YHe = float(_bbn.Y_He(pars.ombh2, 6.14*DELTA))
        if DELTA or self.beta is not None:
            a, w = table(pars, self.beta)
            model = "ppf" if self.beta is not None else "fluid"
            pars.set_dark_energy_w_a(a, w, dark_energy_model=model)
            if model == "fluid": pars.DarkEnergy.cs2 = CS2
        return pars
class CAMBGF(_Base): beta = None
class CAMBLawGF(_Base): beta = 0.5
