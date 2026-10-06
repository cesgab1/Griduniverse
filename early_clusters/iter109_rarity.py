import numpy as np, camb
from scipy.integrate import quad, simpson
h = 0.6736; Om = 0.3153
pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, As=2.1e-9, ns=0.9649)
pars.set_matter_power(redshifts=[0.0], kmax=50.0); res = camb.get_results(pars)
kh, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=50, npoints=600); pk = pk[0]
rho_m = 2.775e11*Om                      # h^2 Msun/Mpc^3 -> in (Msun/h)/(Mpc/h)^3
def sigma(M):                             # M in Msun/h, z = 0
    R = (3*M/(4*np.pi*rho_m))**(1/3); x = kh*R
    W = 3*(np.sin(x) - x*np.cos(x))/x**3
    return np.sqrt(simpson(pk*W**2*kh**2, x=kh)/(2*np.pi**2))
def D(z):                                  # growth factor normalised to 1 today (LCDM)
    a = 1/(1+z); f = lambda a: 1/(a*np.sqrt(Om/a**3 + 1 - Om))**3
    g = lambda a: 2.5*Om*np.sqrt(Om/a**3 + 1 - Om)*quad(f, 1e-6, a)[0]
    return g(a)/g(1)
def dndlnM(M, z):                          # Sheth-Tormen, (h/Mpc)^3
    s = sigma(M)*D(z); e = 1e-3; ds = (np.log(sigma(M*(1+e))*D(z)) - np.log(s))/np.log(1+e)
    nu = 1.686/s; A, a, p = 0.3222, 0.707, 0.3
    fnu = A*np.sqrt(2*a/np.pi)*(1 + (a*nu**2)**-p)*nu*np.exp(-a*nu**2/2)
    return rho_m/M*fnu*abs(ds)
def N_above(Mmin, z1, z2, area_arcmin2):
    area = area_arcmin2/(60**2)*(np.pi/180)**2   # sr
    def dVdz(z): return res.comoving_radial_distance(z)**2*h**3*(2.998e5/(res.hubble_parameter(z)))*h  # (Mpc/h)^3 per sr per dz
    lnM = np.linspace(np.log(Mmin), np.log(Mmin*100), 40)
    def n(z): return simpson([dndlnM(np.exp(l), z) for l in lnM], x=lnM)
    zz = np.linspace(z1, z2, 9)
    return area*simpson([n(z)*dVdz(z) for z in zz], x=zz)
out = []
for Msun in (2e13, 1e13, 5e12):
    N = N_above(Msun*h, 5.2, 6.2, 175)
    out.append(f"halos >= {Msun:.0e} Msun at 5.2<z<6.2 in ~175 arcmin^2: expected {N:.1e}")
out.append(f"sigma(2e13 Msun) at z=5.7: {sigma(2e13*h)*D(5.7):.3f} -> peak height nu = {1.686/(sigma(2e13*h)*D(5.7)):.1f}")
txt = "\n".join(out); print(txt); open("iter109_rarity.txt", "w").write(txt + "\n")
