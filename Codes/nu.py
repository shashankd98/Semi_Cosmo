import numpy as np 
from matplotlib import pyplot as plt  
from scipy.integrate import quad
from scipy.optimize import newton

h=0.7
H0=70.e5/3.086e24
Om0 =0.3089
Ol0=1.-Om0
O0=Ol0+Om0
Ob=0.0469
Msun=2.e33
Gamma=h*Om0*np.exp(-Ob*(1+np.sqrt(2*h)/Om0))
sigma8=0.82

def Ez(z):
    return np.sqrt(Ol0+Om0*(1+z**3))

def Om(z):
    return Om0*(1+z)**3/Ez(z)**2

def Ol(z):
    return Ol0/Ez(z)**2

def g(z):
    return 5./2.*Om(z)*(Om(z)**(4./7.)-Ol(z)+(1.+Om(z)/2.)*(1.+Ol(z)/70.))**(-1)

# def D(z):
#     return g(z)/(1.+z)

def D(z):
    return H(z)/H0*(quad(integrand,z,np.inf)[0]/quad(integrand,0,np.inf)[0])

# def D_norm(z):
#     return D(z)/D(0)

def del0c(z):
    omega=Om0*(1+z)**3/(Ol0+Om0*(1+z)**3)
    p=0.0055
    return 0.15*(12*np.pi)**(2/3)*omega**p

def delc(z):
    return del0c(z)/D(z)

def f(u):
    return 64.087*(1 + 1.074*u**0.3 - 1.581*u**0.4 + 0.954*u**0.5 - 0.185*u**0.6)**(-10)

def sigma(M):
    u8=32*Gamma
    u=3.804e-4*Gamma*(M/(Msun/h*Om0))**(1/3)
    return sigma8*f(u)/f(u8)

def rhs(M0):
    factor=0.254
    return delc(0)+0.477*np.sqrt(2*(sigma(factor*M0)**2-sigma(M0)**2))

def solver(zf,M0):
    return delc(zf)-rhs(M0)

def z_f(M0):

    return newton(solver,1,args=[M0])

def M(z,M0):
    zf=z_f(M0)
    # print(zf)
    nu=1.211+1.858*np.log10(1+zf)+0.308*Ol0**2-0.032*np.log10(M0/(1e11/h*Msun))
    rhs=-0.301*(np.log10(1+z)/np.log10(1+zf))**nu
    # print(z,rhs)
    return M0*10**rhs

def H(z):
    return H0*Ez(z)

def integrand(z):
    return (1+z)/H(z)**3

def D1(z):
    return H(z)/H0*(quad(integrand,z,np.inf)[0]/quad(integrand,0,np.inf)[0])
    # return 5*Om0/2*H0**2*H(z)*quad(integrand,z,np.inf)[0]


def main():
    zs=[0,0.5,1,2,4,6]
    Mass=1e10*Msun*np.logspace(0,5.5,1000)
    nu=np.zeros(1000)
    for z in zs:
        for i in range(1000):
            sig=sigma(Mass[i])
            nu[i]=1.686/sig/D(z)
        plt.semilogy(nu,Mass/Msun,label="z="+str(z))
    plt.xlim(0,4)
    plt.legend()    
    plt.grid()
    plt.show()

    # z=np.linspace(0,10,100)
    # fac=np.zeros(100)
    # for i in range(100):
    #     fac[i]=D(z[i])
    # plt.plot(z,fac)
    # plt.show()

if __name__ == "__main__":
    main()

