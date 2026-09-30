"""Full grid model: dark-energy law (beta, env GRID_BETA) + electron heavier at recombination (VC_ME) relaxing at VC_ZS."""
import os, ctypes
from camb.baseconfig import camblib
from law_camb import CAMBLaw
camblib.set_vconst.argtypes = [ctypes.c_double]*2; camblib.set_vswitch.argtypes = [ctypes.c_double]*2
ME = float(os.environ.get("VC_ME", 1)); ZS = float(os.environ.get("VC_ZS", -1))
class CAMBGrid(CAMBLaw):
    beta: float = float(os.environ.get("GRID_BETA", 0.5))
    def set(self, params_values_dict, state):
        camblib.set_vconst(1.0, ME); camblib.set_vswitch(ZS, 30.0)
        pars = super().set(params_values_dict, state)
        if pars is not None: pars.Recomb.use_rosenbrock = False
        return pars
