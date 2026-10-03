"""Profile chi2 vs link-stretch exponent s; best s and 1-sigma range per SN set (parabola through minimum region)."""
import json, numpy as np, collections
R = collections.defaultdict(dict); L = {}
for line in open("scan_raw.txt"):
    sn, q, tag, _, js = line.split(" ", 4); d = json.loads(js)
    if tag == "LCDM": L[(sn, q)] = d["LCDM"]["chi2"]
    else: R[(sn, q)][float(tag[2:])] = (d["LAW"]["chi2"], d["LAW"]["params"])
out = ["Link-stretch exponent s: physical link length ∝ a^s; law d ln rho_DE/d ln a = (q + 1 - s)/2",
       "s = 1: links stretch with space (Claim 1, rho ∝ adot^-1/2);  s = 0: cells added, fixed link length (rho ∝ H^-1/2)",
       "Fixes Oct 2026 (independent check): sound horizon now integrated from exactly z*, l_A recalibrated to CAMB; exact-q form now\nself-consistent for any s. Significance of s = 0 read directly from the scan values (not extrapolated).",
       "Data: DESI DR2 BAO + Planck 2018 distance priors + one supernova set. Delta chi2 relative to LCDM (same data). No free DE parameter at fixed s.", ""]
for key in sorted(R):
    ss = np.array(sorted(R[key])); c = np.array([R[key][x][0] for x in ss]) - L[key]
    out.append(f"{key[0]:9s} q-form {key[1]:6s}: " + "  ".join(f"s={x:g}:{y:+.1f}" for x, y in zip(ss, c)))
    i = np.argmin(c); sel = slice(max(i - 2, 0), i + 3); p = np.polyfit(ss[sel], c[sel], 2)
    sb = -p[1]/(2*p[0]); sig = 1/np.sqrt(p[0]); cmin = np.polyval(p, sb)
    out.append(f"           best s = {sb:.2f} ± {sig:.2f}  (Delta chi2 vs LCDM at best {cmin:+.1f});  s=1 is {abs(1-sb)/sig:.1f} sigma from best; s=0: sqrt(chi2(s=0)-min) = {np.sqrt(c[ss == 0][0] - cmin):.1f} sigma")
txt = "\n".join(out); print(txt); open("stretch_scan.txt", "w").write(txt + "\n")
