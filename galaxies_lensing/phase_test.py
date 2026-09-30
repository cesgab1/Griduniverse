"""Phase-change dark fluid (superfluid dark matter, Berezhiani & Khoury 2015) — necessary condition with our data.
Coherent (superfluid) where the particles' quantum wavelengths overlap:  n λ³ > 2.612, n = ρ/m, λ = h/(√(2π) m σ)
  ->  coherent if  m < m_crit(ρ, σ) = [ ρ h³ / (2.612 (2π)^{3/2} σ³) ]^{1/4}
Galaxies (where the slack rule works) must be coherent: m < min m_crit(galaxies)   -> m_max
Clusters (where the 3x extra pull appears) must be melted: m > max m_crit(clusters) -> m_min
One particle mass works only if m_min < m_max.  ρ = dark (missing) density inferred from the data; σ = v_c/√2."""
import numpy as np, glob
G=6.674e-11; Msun=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; h=6.626e-34; eV=1.783e-36  # kg per eV/c^2
def mcrit(rho, sig):   # rho kg/m^3, sig m/s -> eV
    return (rho*h**3/(2.612*(2*np.pi)**1.5*sig**3))**0.25/eV
out={}
# --- 1. SPARC outer radii
Yd,Yb=0.5,0.7; m_sp=[]; rows=[]
for f in glob.glob("/home/claude/sparc/r1/*_rotmod.dat"):
    d=np.loadtxt(f); d=d[None] if d.ndim==1 else d
    r,V,eV_,Vg,Vd,Vb=d.T[:6]; ok=(V>0)&(eV_/V<0.1)
    if ok.sum()<4: continue
    r,V,Vg,Vd,Vb=r[ok],V[ok],Vg[ok],Vd[ok],Vb[ok]
    Vbar2=Vg*np.abs(Vg)+Yd*Vd**2+Yb*Vb**2
    Mmiss=np.clip(V**2-Vbar2,0,None)*1e6*r*kpc/G
    dM=np.gradient(Mmiss, r*kpc); rho=dM/(4*np.pi*(r*kpc)**2)
    i=-2; sig=V[-1]*1e3/np.sqrt(2)
    if rho[i]>0: m_sp.append(mcrit(rho[i],sig)); rows.append((f.split('/')[-1][:-11],r[i],V[-1],rho[i],m_sp[-1]))
m_sp=np.array(m_sp)
print(f"SPARC ({len(m_sp)} galaxies, outermost radii ~{np.median([x[1] for x in rows]):.0f} kpc): m_crit median {np.median(m_sp):.2f} eV, 5th percentile {np.percentile(m_sp,5):.2f} eV, min {m_sp.min():.2f} eV")
# --- 2. KiDS lensing rotation curves (isolated galaxies), out to 1 and 3 Mpc
G4=4.52e-30
print("KiDS isolated galaxies (lensing rotation curves):")
m_k1=[]
for k in range(1,5):
    d=np.loadtxt(f"/home/claude/kids/Fig-3_Lensing-rotation-curves_Massbin-{k}.txt"); R=d[:,0]; esd=d[:,1]/d[:,4]
    v=np.sqrt(np.clip(4*G4*esd*R*1e6,0,None))*3.086e13     # km/s (B21 eq. 23)
    rv2=R*Mpc*(v*1e3)**2; rho=np.gradient(rv2,R*Mpc)/(4*np.pi*G*(R*Mpc)**2)
    for Rt in (0.3,1.0,3.0):
        j=np.argmin(np.abs(R-Rt)); rr=np.interp(Rt,R,np.clip(rho,1e-40,None)); vv=np.interp(Rt,R,v)
        mc=mcrit(rr,vv*1e3/np.sqrt(2))
        if Rt<=1.0: m_k1.append(mc)
        print(f"   mass bin {k}: R={Rt:.1f} Mpc  v≈{vv:.0f} km/s  ρ≈{rr*1e-3:.1e} g/cm³  m_crit {mc:.2f} eV")
# --- 3. clusters (CLASH stacked NFW, Umetsu+2016-like: M200c 1.2e15 Msun, c 3.7) and a range
rhoc=3*(70e3/Mpc)**2/(8*np.pi*G)
def nfw(M200,c,r):
    r200=(3*M200*Msun/(4*np.pi*200*rhoc))**(1/3); rs=r200/c
    dc=200/3*c**3/(np.log(1+c)-c/(1+c)); x=r/rs
    rho=rhoc*dc/(x*(1+x)**2); M=4*np.pi*rhoc*dc*rs**3*(np.log(1+x)-x/(1+x)); v=np.sqrt(G*M/r)
    return rho, v, r200
print("Clusters (NFW, where the cluster relation departs from the slack rule, 0.2-1.5 Mpc):")
m_cl=[]
for M200,c in [(5e14,4.0),(1.2e15,3.7),(2e15,3.5)]:
    for rM in (0.2,0.5,1.0,1.5):
        rho,v,r200=nfw(M200,c,rM*Mpc); mc=mcrit(rho,v/np.sqrt(2)); m_cl.append(mc)
        print(f"   M200 {M200:.0e}: r={rM:.1f} Mpc  σ≈{v/np.sqrt(2)/1e3:.0f} km/s  ρ≈{rho*1e-3:.1e} g/cm³  m_crit {mc:.2f} eV")
# --- 4. groups (need coherent at r <= 0.2 Mpc per KiDS inner-radius result)
print("Groups (NFW M200 1e13-1e14, r = 0.1-0.2 Mpc, must stay coherent):")
m_gr=[]
for M200,c in [(1e13,7),(3e13,6),(1e14,5)]:
    for rM in (0.1,0.2):
        rho,v,_=nfw(M200,c,rM*Mpc); mc=mcrit(rho,v/np.sqrt(2)); m_gr.append(mc)
        print(f"   M200 {M200:.0e}: r={rM:.1f} Mpc σ≈{v/np.sqrt(2)/1e3:.0f} km/s  m_crit {mc:.2f} eV")
m_max_gal=min(np.percentile(m_sp,5), min(m_k1)); m_max_grp=min(m_gr); m_min=max(m_cl)
print(f"\nNeed m < {m_max_gal:.2f} eV (galaxies, incl. KiDS to 1 Mpc), m < {m_max_grp:.2f} eV (groups), and m > {m_min:.2f} eV (clusters, all radii 0.2-1.5 Mpc)")
print("window EXISTS" if m_min<min(m_max_gal,m_max_grp) else "NO window")
# restricted: clusters only need melting where the deficit is largest (outskirts, r>=0.5 Mpc)

# ---- relaxed reading: slack verified around galaxies to 0.3 Mpc (KiDS inner, groups to 0.2 Mpc); cluster deficit strongest in outskirts (>=1 Mpc)
m_max=min(np.percentile(m_sp,5), 3.02, m_max_grp); m_min_out=max(mcrit(*[x for x in (nfw(M,c,1.0*Mpc)[0],)], nfw(M,c,1.0*Mpc)[1]/np.sqrt(2)) for M,c in [(5e14,4.0),(1.2e15,3.7),(2e15,3.5)])
print(f"\nRelaxed: coherent galaxies to 0.3 Mpc & groups to 0.2 Mpc -> m < {m_max:.2f} eV; melted cluster outskirts (>=1 Mpc) -> m > {m_min_out:.2f} eV")
print("window:", f"{m_min_out:.2f}-{m_max:.2f} eV" if m_min_out<m_max else "none")
# where does each system switch phase, for masses inside the window?
def switch_radius(profile, m):
    rs=np.geomspace(0.01,3,400)
    mc=np.array([profile(r) for r in rs])
    melted=mc<m
    return rs[np.argmax(melted)] if melted.any() else np.inf
for m in (2.0,2.5,3.0):
    print(f"\nm = {m} eV: phase-switch radius (inside coherent, outside melted)")
    for M200,c in [(1e13,7),(1e14,5),(5e14,4),(1.2e15,3.7)]:
        r=switch_radius(lambda rr: mcrit(nfw(M200,c,rr*Mpc)[0], nfw(M200,c,rr*Mpc)[1]/np.sqrt(2)), m)
        print(f"   halo M200 {M200:.0e} (σ≈{nfw(M200,c,0.3*Mpc)[1]/np.sqrt(2)/1e3:.0f} km/s): switch at {r:.2f} Mpc")
    for k in range(1,5):
        d=np.loadtxt(f"/home/claude/kids/Fig-3_Lensing-rotation-curves_Massbin-{k}.txt"); R=d[:,0]; esd=d[:,1]/d[:,4]
        v=np.sqrt(np.clip(4*G4*esd*R*1e6,0,None))*3.086e13; rv2=R*Mpc*(v*1e3)**2
        rho=np.clip(np.gradient(rv2,R*Mpc)/(4*np.pi*G*(R*Mpc)**2),1e-40,None)
        mc=mcrit(rho,v*1e3/np.sqrt(2)); melted=mc<m
        print(f"   KiDS isolated galaxies bin {k}: switch at {R[np.argmax(melted)] if melted.any() else np.inf:.2f} Mpc")
