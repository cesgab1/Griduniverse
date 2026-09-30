import json,numpy as np
P='/home/claude/fullfit/pact/'
rows=[('ΛCDM','lcdm'),('β=0 (constant tension) + electron 1.12%','b0E'),('β=1 (coherent) + electron 0.17%','b1E'),('grid, electron from energy budget +0.43% (β law, self-consistent)','gridE043'),('grid, electron from energy budget +1.12% (Λ budget)','gridE112'),('heavier electron only','me'),('dark-energy law β=½ only','b05'),('full grid (β=½ + electron)','grid')]
full=json.load(open(P+'vc_me1.00.json'))   # ΛCDM fitted to everything, for reference
print("reference ΛCDM fitted to all data: BAO χ²",round(full['parts'].get('chi2__bao.desi_dr2',0),1),"SN χ²",round(full['parts'].get('chi2__sn.desdovekie',0),1))
out=[]
for lab,t in rows:
    f=json.load(open(P+f'fwd_{t}.json')); e=json.load(open(P+f'ev_{t}.eval.json'))
    p=f['params']; d=f['derived']; h=p['H0']/100; Om=(p['ombh2']+p['omch2']+0.00064)/h**2
    S8=d['sigma8']*np.sqrt(Om/0.3)
    r=dict(model=lab,cmb=f['chi2'],H0=p['H0'],Om=Om,S8=S8,age=d['age'],bao=e['chi2__bao.desi_dr2'],sn=e['chi2__sn.desdovekie'])
    out.append(r)
b=out[0]
for r in out:
    print(f"{r['model']:30s} CMBχ² {r['cmb']-b['cmb']:+6.1f} | H0 {r['H0']:5.2f} Ωm {r['Om']:.3f} S8 {r['S8']:.3f} age {r['age']:.2f} Gyr | predicted DESI BAO χ² {r['bao']:6.1f} ({r['bao']-b['bao']:+5.1f}) DES SN χ² {r['sn']:7.1f} ({r['sn']-b['sn']:+5.1f})")
json.dump(out,open('/tmp/claude-0/-home-claude-griduniverse/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/forward_summary.json','w'),indent=1)
