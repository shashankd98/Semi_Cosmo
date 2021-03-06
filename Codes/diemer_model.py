#Playing aroun with Diemer-Kravtsov 2014 DM profiles & associated potential
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.special import gammainc
from scipy.special import gamma
def Ez(z):
    return np.sqrt( Om*(1.+z)**3+(1.-Om) )
def Ezsqr(z):
    return Om*(1.+z)**3+(1.-Om)
def delc(x): #x=Om(z)-1
    return 18.*np.pi**2 + 82.*x - 39.*x*x
def d_Ein(r):
    return ds*np.exp( -2.*((r/rs)**alp-1.)/alp )
def d_NFW(r):
    c=R200c/rs
    return 200.*c*c*c*dcr/(3.*(np.log(1.+c)-c/(1.+c)))/(r/rs*(1+r/rs)**2)
def f_tr(r):
    return (1. + (r/rt)**bet)**(-gam/bet)
def d_out(r):
    return dm*(be*(r/(5.*R200m))**-se + 1.)
def d_rho(r):
    #return d_Ein(r)
    #return d_NFW(r)
    #return d_Ein(r)*f_tr(r)+d_out(r)
    return d_NFW(r)*f_tr(r)+d_out(r)
def Menc_NFW(r):
    x=r/R200c; c=R200c/rs
    return M200*x*( -c/(1.+c*x) + np.log(1.+c*x)/x )/( np.log(1.+c) - c/(1.+c) )
def Menc_Ein(r):
    return 2.*np.pi*ds*rs*rs*rs*np.exp(2./alp)*(alp/2.)**(-1.+3./alp)*gamma(3./alp)*gammainc(3./alp,2.*(r/rs)**alp/alp)
def Phi_NFW(r):
    x=r/R200c; c=R200c/rs
    return -G*M200*np.log(1.+c*x)/(r*(np.log(1.+c)-c/(1.+c)))
def dMbydr(M,r):
    return 4.*np.pi*r*r*d_rho(r)
def Menc(r):
    #sol = odeint( dMbydr, Menc_Ein(1.e-2*R200m), r )
    sol = odeint( dMbydr, Menc_NFW(1.e-2*R200m), r )
    return sol[:,0]
def Menc_phi(r):
    rarray = np.array([1.e-2*R200m, r])
    sol = odeint( dMbydr, Menc_NFW(1.e-2*R200m), rarray )
    #sol = odeint( dMbydr, Menc_Ein(1.e-2*R200m), r )
    return sol[1,0]
def dPhibydr(Phi,r):
    return G*Menc_phi(r)/(r*r)
def sigmoid(x):
    return 1 / (1 + np.exp(-x))
global Om, ds, alp, rs, bet, rt, gam, dm, se, dcr, R200m, R200c, M200, G
#cosmological, astrophysical parameters
H0 = 70.e5/3.086e24; G=6.67e-8; z=0.0; Om = 0.3089; dcr0 = 3.0*H0*H0/(8.*np.pi*G)
dcr=dcr0*Ezsqr(z); Msun=2.e33; dm=dcr0*Om*(1.+z)**3
#halo parameters; R200c used for NFW, R200m for Ein
M200 = 1.e12*Msun; 
#M200 = 1.e12*Msun
c=4.0; nu=4.0; #
#c=4.0; nu=0.5 #parameters for low mass halos
alp=0.155+0.0095*nu*nu 
R200c = (3.*M200/(800.*np.pi*dcr))**(1./3.); R200m = (3.*M200/(800.*np.pi*dm))**(1./3.)
rs=R200c/c; R200m=(Ezsqr(z)/(Om*(1.+z)**3))**(1./3.)*R200c
ds=M200/(2.*np.pi*rs*rs*rs*np.exp(2./alp)*(alp/2.)**(3./alp-1.)*gamma(3./alp)*gammainc(3./alp,2.*(R200m/rs)**alp/alp))
#parameters for outer profile
bet=6.0; gam=4.0; #bet=4.0; gam=8.0;
be = 1.5; se=1.5; rt=(1.9-0.18*nu)*R200m 
r=R200m*np.logspace(-2,3,500); 
#sol = odeint( dMbydr, Menc_Ein(1.e-2*R200m), r )
#sol = odeint( dMbydr, Menc_NFW(1.e-2*R200m), r )
sol_phi = odeint( dPhibydr, Phi_NFW(1.e-2*R200m), r)
g = G*Menc(r)/(r*r); 
tff = np.sqrt(2.*r/g)
#plt.plot(r/R200c, -Phi_NFW(r), r/R200c, -sol_phi[:,0], r/R200c, sol_phi[:,0])
#plt.plot(r/R200m,sol[:,0])
#plt.plot(r/R200c, g,'-',r/R200c, G*Menc_NFW(r)/(r*r), r[:-1]/R200c, np.diff(sol_phi[:,0])/np.diff(r), r[:-1]/R200c, np.diff(Phi_NFW(r))/np.diff(r))
#plt.plot(r/R200m, g,'-',r/R200m, G*Menc_Ein(r)/(r*r))
#plt.plot(r/R200m, tff/np.pi/1.e16,'o')
rscale1=10*R200c; r1=10*R200c; p1=1./3.; rscale2=10.*R200c; r2=10.*R200c; p2=1./3.

f1 = Phi_NFW(r); f2 = (4.*np.pi*G)*d_rho(R200c*100.)*(be*(5.*R200m)**(se)*(r**(2.-se)/((3.-se)*(2.-se)))+r**2/6.)

w1 =1.-sigmoid(1.2*np.log(r/rt)); w2 =1.-w1
phi_fit = f1*w1+f2*w2
#plt.plot(r/R200c,w1,r/R200c,w2,r/R200c,w1+w2)

# plt.subplot(311)
# plt.plot(r/R200c, sol_phi[:,0], r/R200c, -sol_phi[:,0])
# plt.plot(r/R200c, phi_fit, r/R200c, -phi_fit)
# plt.plot(r/R200c, -Phi_NFW(r))
# plt.ylabel(r'$\Phi (cm^2s^{-2})$')
# plt.grid()
# #plt.axis([1, 100, 1e13, 1e16])
# plt.xscale('log'); plt.yscale('log')
# plt.show()
# plt.subplot(312)
plt.plot(r,phi_fit)
plt.show()

g_fit = np.diff(phi_fit)/np.diff(r)
'''
plt.plot(r/R200c,g,r[:-1]/R200c,g_fit, r/R200c, G*Menc_NFW(r)/r**2)
plt.grid()
plt.xscale('log'); plt.yscale('log')
#plt.axis([1, 100, 1e-13, 1e-8])
plt.ylabel(r'$g(cm s^{-2})$')
plt.show()
'''

plt.subplot(211)
rm1=r[:-1]; rm2=rm1[:-1]
rho_fit = np.diff(rm1*rm1*g_fit)/np.diff(rm1*rm1*rm1/3.)/(4.*np.pi*G)
plt.plot(r,d_rho(r))
plt.plot(rm1[:-1],rho_fit)
#plt.plot(r[:], sol_phi)
#plt.plot(r[:],phi_fit)
plt.grid()
plt.xscale('log'); plt.yscale('log')
#plt.axis([1, 100, 1e-30, 1e-27])
plt.ylabel(r'$\rho (g cm^{-3})$')
plt.xlabel('$r/R_{200c}$')
#plt.show()
plt.subplot(212)
plt.plot(r[:-1]/R200c,np.diff(np.log(d_rho(r)))/(np.diff(np.log(r))))
plt.plot(rm2[:-1]/R200c,np.diff(np.log(rho_fit))/(np.diff(np.log(rm2))))
plt.grid()
plt.xscale('log')
plt.show()
