"""T2 sensitivity: how the verdict on each named alpha depends on the frequency tolerance (free choice)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, qnm
mode = qnm.modes_cache(s=-2, l=2, m=2, n=0)
chis = np.linspace(0.67, 0.69, 41); Fs = np.array([mode(a=x)[0].real for x in chis])
named = {"4ln2 (one bit)":4*np.log(2), "4ln3 (Hod)":4*np.log(3), "8ln2":8*np.log(2), "8pi":8*np.pi}
alphas = np.round(np.arange(0.5, 120.0001, 0.01), 2)
lines = []
for tol in (0.020, 0.024, 0.029, 0.040, 0.050):
    ok = np.zeros(len(alphas), bool)
    for ch, F in zip(chis, Fs):
        s = np.sqrt(1-ch**2); step = alphas*s/(16*np.pi*(1+s)); x = (F-ch/(1+s))/step
        ok |= np.abs(x-np.round(x))*step <= tol*F
    v = "  ".join(f"{k}:{'ok' if ok[np.argmin(abs(alphas-a))] else 'X'}" for k,a in named.items())
    lines.append(f"tol {tol*100:4.1f}%: smallest excluded alpha {alphas[~ok].min():6.2f}; excluded fraction {(~ok).mean()*100:3.0f}%  | {v}")
print("\n".join(lines)); open("RESULT_123_t2_sensitivity.txt","w").write("GW250114, spin 0.67-0.69\n"+"\n".join(lines)+"\n")
