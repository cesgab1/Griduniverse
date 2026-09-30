"""CAMB (patched: G(a)/G0 = 1 + delta/(1+(a/a_t)^4) in every Einstein equation) for Cobaya.
Helium from BBN with the boost included (equivalent dN_eff = 6.14 delta at neutron freeze-out)."""
import os, ctypes
import camb.bbn as bbn
from camb.baseconfig import camblib
from cobaya.theories.camb import CAMB
from law_camb import CAMBLaw
_set = camblib.set_gboost; _set.argtypes = [ctypes.c_double, ctypes.c_double]
DELTA = float(os.environ.get("GDELTA", 0)); AT = float(os.environ.get("GAT", 1e-4))
_bbn = bbn.get_predictor()
def _fix(pars):
    _set(DELTA, AT)
    if pars is not None:
        pars.YHe = float(_bbn.Y_He(pars.ombh2, 6.14 * DELTA))
    return pars
class CAMBG(CAMB):
    def set(self, params_values_dict, state):
        return _fix(super().set(params_values_dict, state))
class CAMBLawG(CAMBLaw):
    def set(self, params_values_dict, state):
        return _fix(super().set(params_values_dict, state))
