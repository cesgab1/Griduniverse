"""Iteration 77: are neutrinos dark energy's clock? See PREREG_77.md."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid as ctz
dm21, dm31 = 7.49e-5, 2.513e-3            # eV^2
Tnu0 = 1.945*8.617e-5                      # eV
pbar0 = 3.15*Tnu0
ODE, Ok, Or = 0.685, 0.0023, 9.1e-5; Om = 1 - ODE - Ok - Or; tH = 977.8/67.4
x = np.linspace(-9, 0, 20001)
def E1(xi):
    z1 = np.exp(-xi); base = Om*z1**3 + Or*z1**4 + Ok*z1**2
    return np.exp(brentq(lambda l: np.exp(2*l) - base - ODE*(np.exp(l)*np.exp(xi))**-0.5, -80, 80))
E = np.array([E1(v) for v in x]); t = (ctz(1/E, x, initial=0) + 1/(2*E[0]))*tH
t_of_z = lambda z: np.interp(-np.log(1+z), x, t)
out = ["Iteration 77: neutrinos as dark energy's clock (measured inputs). PREREG_77.md", ""]
m2min, m3min = np.sqrt(dm21), np.sqrt(dm31)
out.append("N1 when the measured-scale neutrinos slow down (criterion 1 | criterion 2):")
for name, m in (("m2 >= 8.65 meV", m2min), ("m3 >= 50.1 meV", m3min)):
    z1, z2 = m/pbar0 - 1, np.sqrt(3)*m/pbar0 - 1   # careful: criterion 2 needs p = sqrt(3) m -> earlier? p falls with time, so p = sqrt3 m happens EARLIER
    out.append(f"   {name}: z = {z1:.1f} (t = {t_of_z(z1)*1e3:.0f} Myr) | z = {z2:.1f} (t = {t_of_z(z2)*1e3:.0f} Myr)  vs onset 7.47 Gyr -> MISS")
out.append("")
out.append("N2 lightest neutrino mass needed to slow down 'now' (a FIT, not a prediction):")
for label, z in (("acceleration onset z=0.67", 0.67), ("dark energy passes matter z=0.32", 0.32)):
    for crit, f in (("crit 1", 1.0), ("crit 2", 1/np.sqrt(3))):
        m1 = pbar0*(1+z)*f
        s_no = m1 + np.sqrt(m1**2 + dm21) + np.sqrt(m1**2 + dm31)
        m3 = m1; s_io = m3 + np.sqrt(m3**2 + dm31) + np.sqrt(m3**2 + dm31 + dm21)
        out.append(f"   {label}, {crit}: m1 = {m1*1e3:.2f} meV -> sum {s_no*1e3:.1f} meV (normal; bound 64.2 -> "
                   f"{'ok' if s_no < 0.0642 else 'EXCLUDED'}) | inverted {s_io*1e3:.1f} meV -> {'ok' if s_io < 0.0642 else 'EXCLUDED'}")
s0 = np.sqrt(dm21) + np.sqrt(dm31)
out.append(f"   for comparison m1 = 0: sum {s0*1e3:.1f} meV -> the 'now' prediction differs from zero by < 1 meV in the sum (unmeasurable today)")
out.append("")
de = 2.24e-3
out.append(f"N3 energy scales: dark energy {de*1e3:.2f} meV vs 8.65 meV (x{m2min/de:.1f}) and 50.1 meV (x{m3min/de:.1f}) -> "
           + ("PASS" if min(m2min/de, m3min/de) < 2 else "FAIL"))
h = 0.674; Onu = 0.0595/93.14/h**2
out.append(f"N4 neutrino energy today (sum 59.5 meV): {Onu*100:.2f}% of the total; dark energy is {ODE/Onu:.0f} x larger -> neutrinos cannot BE dark energy.")
txt = "\n".join(out); print(txt); open("iter77_neutrino_clock.txt", "w").write(txt + "\n")
