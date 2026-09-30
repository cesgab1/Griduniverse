"""Minimal cobaya Theory serving a tabulated background (H(z), D_A(z), rdrag) for BAO + SN likelihoods."""
import numpy as np
from cobaya.theory import Theory
C = 299792.458
class BG(Theory):
    params = {"dummy": None, "rdrag": {"derived": True}}
    state_bg = None   # set externally: dict(z=array ascending, H=km/s/Mpc array, rdrag=Mpc)
    def get_can_provide(self): return ["Hubble", "angular_diameter_distance", "comoving_radial_distance"]
    def get_can_provide_params(self): return ["rdrag"]
    def must_provide(self, **req): return {}
    def calculate(self, state, want_derived=True, **params):
        b = BG.state_bg; z = b["z"]; H = b["H"]
        chi = np.concatenate([[0], np.cumsum(0.5*(C/H[1:]+C/H[:-1])*np.diff(z))])
        state["z"] = z; state["H"] = H; state["chi"] = chi; state["derived"] = {"rdrag": b["rdrag"]}
    def get_Hubble(self, z, units="km/s/Mpc"):
        H = np.interp(z, self.current_state["z"], self.current_state["H"])
        return H if units == "km/s/Mpc" else H/C
    def get_comoving_radial_distance(self, z): return np.interp(z, self.current_state["z"], self.current_state["chi"])
    def get_angular_diameter_distance(self, z): return self.get_comoving_radial_distance(z)/(1+np.asarray(z))
