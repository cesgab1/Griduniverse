# Third requirement: the same mix inside galaxies. If a fraction f of the dark ingredient gathers into structures, galaxies get a halo too.
# Halo pull estimated from the data themselves: in the standard picture the halo supplies g_obs - g_bar; with the mix it supplies f*(g_obs - g_bar),
# and slack then acts on the total. Compared with what is measured (KiDS isolated lensing; SPARC rotation curves).
import numpy as np, glob
c=2.998e8; Mpc=3.0857e22; H0=67.73e3/Mpc; a0=c*H0/6
nu = lambda y: 1/(-np.expm1(-np.sqrt(y)))
G4=4.52e-30; pcm=3.086e16; D="/home/claude/kids/"
d=np.loadtxt(D+"Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt"); gbK=d[:,0]; goK=4*G4*d[:,1]/d[:,4]*pcm
Yd,Yb,kpc=0.5,0.7,3.0857e19; GB=[];GO=[]
for f in glob.glob("/home/claude/sparc/r1/*_rotmod.dat"):
    x=np.loadtxt(f); x=x[None] if x.ndim==1 else x
    r,V,eV,Vg,Vd,Vb=x.T[:6]; gb=(Vg*np.abs(Vg)+Yd*Vd**2+Yb*Vb**2)*1e6/(r*kpc); go=V**2*1e6/(r*kpc)
    ok=(V>0)&(eV/V<0.1)&(gb>0); GB+=list(gb[ok]); GO+=list(go[ok])
GB=np.array(GB); GO=np.array(GO)
for lab,gb,go in [("SPARC rotation curves",GB,GO),("KiDS isolated lensing",gbK,goK)]:
    print(f"\n{lab}: predicted/measured pull (median) by strength of visible pull")
    edges=[1e-13,1e-12,1e-11,1e-10,1e-9]
    for f in [0.0,0.74,0.93]:
        gt=gb+f*np.clip(go-gb,0,None); pred=gt*nu(gt/a0); rat=pred/go
        cells=[]
        for lo,hi in zip(edges[:-1],edges[1:]):
            m=(gb>=lo)&(gb<hi); cells.append(f"{np.median(rat[m]):.2f}" if m.sum()>2 else "  - ")
        print(f"   f = {f:.2f} (slack {'only' if f==0 else '+ halo'}):  g_bar 1e-13..1e-12: {cells[0]}   1e-12..1e-11: {cells[1]}   1e-11..1e-10: {cells[2]}   1e-10..1e-9: {cells[3]}")
