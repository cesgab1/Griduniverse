"""CAMB with the fine-structure constant and electron mass at recombination set by env VC_ALPHA, VC_ME (ratios to today).
Recombination physics scaled (energies ~ alpha^2 m_e, rates, Thomson cross-section ~ alpha^2/m_e^2)."""
import os, ctypes
from camb.baseconfig import camblib
from cobaya.theories.camb import CAMB
camblib.set_vconst.argtypes = [ctypes.c_double, ctypes.c_double]
camblib.set_vswitch.argtypes = [ctypes.c_double, ctypes.c_double]
ZS = float(os.environ.get("VC_ZS", -1)); DZ = float(os.environ.get("VC_DZ", 30))
AL = float(os.environ.get("VC_ALPHA", 1)); ME = float(os.environ.get("VC_ME", 1))
class CAMBVC(CAMB):
    def set(self, params_values_dict, state):
        pars = super().set(params_values_dict, state)
        camblib.set_vconst(AL, ME)
        camblib.set_vswitch(ZS, DZ)
        if pars is not None: pars.Recomb.use_rosenbrock = False
        return pars
