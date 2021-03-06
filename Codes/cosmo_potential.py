#Playing aroun with Diemer-Kravtsov 2014 DM profiles & associated potential
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.special import gammainc
from scipy.special import gamma
from scipy.integrate import quad
from scipy.special import erf

pc=3.086e18

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
def sigmoid(x,ce):
    return (1.+erf(x))/ce
def first_term(r):
    return (d_rho(r)-dm)*r**2 
def second_term(r):
    return (d_rho(r)-dm)*r
def phi_real(r):
    phi=-4*np.pi*G*(quad(first_term,1e-3*R200m,r)[0]/r+quad(second_term,r,9*R200m)[0])
    return phi
def Phi_outer(r):
    phi=(4.*np.pi*G)*dm*(be*(5.*R200m)**(se)*(r**(2.-se)/((3.-se)*(2.-se)))+r**2/6.)
    return phi
global Om, ds, alp, rs, bet, rt, gam, dm, se, dcr, R200m, R200c, M200, G
#cosmological, astrophysical parameters
H0 = 70.e5/3.086e24; G=6.67e-8; z=0.0; Om = 0.3089; dcr0 = 3.0*H0*H0/(8.*np.pi*G)
dcr=dcr0*Ezsqr(z); Msun=2.e33; dm=dcr0*Om*(1.+z)**3
#halo parameters; R200c used for NFW, R200m for Ein
M200 = 1.e14*Msun; 
#M200 = 1.e12*Msun
c=4.0; nu=4.0; #
# c=4.0; nu=0.5 #parameters for low mass halos
alp=0.155+0.0095*nu*nu 
R200c = (3.*M200/(800.*np.pi*dcr))**(1./3.); R200m = (3.*M200/(800.*np.pi*dm))**(1./3.)
rs=R200c/c; R200m=(Ezsqr(z)/(Om*(1.+z)**3))**(1./3.)*R200c

ds=M200/(2.*np.pi*rs*rs*rs*np.exp(2./alp)*(alp/2.)**(3./alp-1.)*gamma(3./alp)*gammainc(3./alp,2.*(R200m/rs)**alp/alp))
#parameters for outer profile
bet=6.0; gam=4.0; #bet=4.0; gam=8.0;
be = 1.5; se=1.5; rt=(1.9-0.18*nu)*R200m 
def main():
    # print(Phi_NFW(R200c)/((4.*np.pi*G)*d_rho(R200c*100.)*(be*(5.*R200m)**(se)*(R200c**(2.-se)/((3.-se)*(2.-se)))+R200c**2/6.)))
    r=1e3*pc*np.logspace(1,5,500); 
    #sol = odeint( dMbydr, Menc_Ein(1.e-2*R200m), r )
    #sol = odeint( dMbydr, Menc_NFW(1.e-2*R200m), r )
    sol_phi = odeint( dPhibydr, Phi_NFW(1.e-2*R200m), r)
    phi_sol=np.zeros(500)
    for i in range(500):
        phi_sol[i]=sol_phi[i]-2./3.*np.pi*G*dm*r[i]*r[i]
        # phi_sol[i]=sol_phi[i]
    g = G*Menc(r)/(r*r); 
    tff = np.sqrt(2.*r/g)
    #plt.plot(r/R200c, -Phi_NFW(r), r/R200c, -sol_phi[:,0], r/R200c, sol_phi[:,0])
    #plt.plot(r/R200m,sol[:,0])
    #plt.plot(r/R200c, g,'-',r/R200c, G*Menc_NFW(r)/(r*r), r[:-1]/R200c, np.diff(sol_phi[:,0])/np.diff(r), r[:-1]/R200c, np.diff(Phi_NFW(r))/np.diff(r))
    #plt.plot(r/R200m, g,'-',r/R200m, G*Menc_Ein(r)/(r*r))
    #plt.plot(r/R200m, tff/np.pi/1.e16,'o')
    rscale1=10*R200c; r1=10*R200c; p1=1./3.; rscale2=10.*R200c; r2=10.*R200c; p2=1./3.

    Omz=Om*pow((1.+z),3.)/Ezsqr(z)
    Phi_cosmo=-H0*H0*Ezsqr(z)*(1.-3.*Omz/2.)*r**2/2.

    g_cosmo=np.diff(-Phi_cosmo)/np.diff(r)

    f1 = Phi_NFW(r); f2 = Phi_outer(r)-2./3.*np.pi*G*dm*r*r
    
    w1 =1-sigmoid(np.log(r/rt),7.); w2 =(1.-w1**4)**(1/8)
    # phi_fit = f1*w1+(f2-2./3.*np.pi*G*dm*r**2)*w2
    phi_fit=f1*w1+f2*w2
    #plt.plot(r/R200c,w1,r/R200c,w2,r/R200c,w1+w2)

    # plt.loglog(r,-Phi_NFW(r))
    # plt.loglog(r,Phi_outer(r))
    # plt.loglog(r,phi_fit)
    # plt.show()

    phi_true=np.zeros(500)
    for i in range(500):
        phi_true[i]=phi_real(r[i])

    # phi_true=phi_true-Phi_cosmo
    # phi_sol=phi_sol-Phi_cosmo

    fig = plt.figure(figsize=(20,10))
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223)
    # ax4=fig.add_subplot(224)

    phi_full=Phi_cosmo+phi_true
    g_full=np.diff(-phi_full)/np.diff(r)


    # ax1.plot(r/R200m,phi_true,label="Only DM")
    # ax1.plot(r/R200m,Phi_cosmo,label='Only Cosmo')
    # ax1.plot(r/R200m,phi_full,label="DM+Cosmo")
    # ax1.set_xlabel("R/$R_{200m}$")
    # ax1.set_ylabel("$\Phi$(R)")
    # ax1.legend()

    # ax1.set_yscale('symlog')
    # ax1.set_xscale('log')

    rm1=r[:-1]; rm2=rm1[:-1]; rm3=rm2[:-1]

    # ax1.loglog(rm1/R200m,g_cosmo)

    g_true=np.diff(-phi_true)/np.diff(r)
    g_fit=np.diff(-phi_fit)/np.diff(r)
    g_sol=np.diff(-phi_sol)/np.diff(r)

    # ax2.plot(rm1/R200m,g_true)
    # ax2.plot(rm1/R200m,g_cosmo)
    # ax2.plot(rm1/R200m,g_full)
    # ax2.set_xlabel("R/$R_{200m}$")
    # ax2.set_ylabel("g(R)")

    # ax2.set_yscale('symlog')
    # ax2.set_xscale('log')

    ax1.loglog(rm1/R200c,-g_true)
    ax1.set_title("DM")
    ax1.set_xlabel("r/R200c")
    ax1.set_ylabel("g (CGS units) (inward)")

    ax2.loglog(rm1/R200c,g_cosmo)
    ax2.set_title("Cosmo")
    ax2.set_xlabel("r/R200c")
    ax2.set_ylabel("g (CGS units)")

    ax3.semilogx(rm1/R200c,g_full)
    ax3.set_title("DM+cosmo")
    # ax3.set_yscale('symlog')
    ax3.set_xlabel("r/R200c")
    ax3.set_ylabel("g (CGS units)")

    # slope=np.diff(np.log(g_true))/np.diff(np.log(rm1))
    # ax3.semilogx(rm2/R200m,slope)
    # ax3.set_xlabel("R/$R_{200m}$")
    # ax3.set_ylabel("$\gamma$")

    # rho_fit = np.diff(rm1*rm1*(-g_fit))/np.diff(rm1*rm1*rm1/3.)/(4.*np.pi*G)
    # rho_true=np.diff(rm1*rm1*(-g_true))/np.diff(rm1*rm1*rm1/3.)/(4.*np.pi*G)
    # rho_sol=np.diff(rm1*rm1*(-g_sol))/np.diff(rm1*rm1*rm1/3.)/(4.*np.pi*G)

    # ax3.loglog(rm2/R200m,rho_fit)
    # ax3.loglog(rm2/R200m,rho_true)
    # ax3.loglog(rm2/R200m,rho_sol)
    # ax3.set_xlabel("R/$R_{200m}$")
    # ax3.set_ylabel(r'$\rho (R)$')

    # slope_fit=np.diff(np.log(rho_fit))/(np.diff(np.log(rm2)))
    # slope_true=np.diff(np.log(rho_true))/(np.diff(np.log(rm2)))
    # slope_sol=np.diff(np.log(rho_sol))/(np.diff(np.log(rm2)))


    # ax4.semilogx(rm3/R200m,slope_fit)
    # ax4.semilogx(rm3/R200m,slope_true)
    # ax4.semilogx(rm3/R200m,slope_sol)
    # ax4.set_xlabel("R/$R_{200m}$")
    # ax4.set_ylabel(r"$\frac{d \log(\rho)}{d \log(r)}$")
    fig.suptitle(r"$M_0=$"+str('%2g'%(M200/Msun)))
    plt.show()
    


if __name__ == "__main__":
    main()
