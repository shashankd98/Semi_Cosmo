import numpy as np 
from matplotlib import pyplot as plt  
from scipy.integrate import quad
from scipy.optimize import newton
import vdb

h=0.7
H0=70.e5/3.086e24
Om0 =0.27
Ol0=1.-Om0
O0=Ol0+Om0
Ob=0.019/h**2
Msun=2.e33
sigma8=1.
Gamma=h*Om0*np.exp(-Ob*(1+np.sqrt(2*h)/Om0))
slope=-0.37865052324944637

def Ez(z):
    return np.sqrt(Ol0+Om0*(1+z**3))

def Om(z):
    return Om0*(1+z)**3/Ez(z)**2

def Ol(z):
    return Ol0/Ez(z)**2

def H(z):
    return H0*Ez(z)

def integrand(z):
    return (1+z)/H(z)**3

def D(z):
    return H(z)*quad(integrand,z,np.inf)[0]

def D_norm(z):
    return D(z)/D(0)

def zf(M0):
    return -0.0064*(np.log10(M0/Msun))**2 + 0.0237*(np.log10(M0/Msun)) + 1.8837

def f(u):
    return 64.087*(1 + 1.074*u**0.3 - 1.581*u**0.4 + 0.954*u**0.5 - 0.185*u**0.6)**(-10)

def q(M0):
    return 4.137*zf(M0)**(-0.9476)

def sigma(M):
    u8=32*Gamma
    u=3.804e-4*Gamma*(M/(Msun/h*Om0))**(1/3)
    return sigma8*f(u)/f(u8)

def fm(M,M0):
    return 1./np.sqrt(sigma(M/q(M0))**2-sigma(M)**2)

def alpha():
    return (1.686*(2./np.pi)**(1./2.)*slope+1.)

def M(z,M0):
    Mz=M0*(1.+z)**(alpha()*fm(M0,M0))*np.exp(-fm(M0,M0)*z)
    return Mz

def main():
    z=np.linspace(0,10,100)
    M0=5*Msun*np.logspace(11,14,4)
    for j in range(4):

        if j==0: color='k'
        elif j==1: color='r'
        elif j==2: color='b'
        elif j==3: color='g'
        elif j==4: color='y'


        yax=np.zeros(100)
        yax2=np.zeros(100)
        for i in range(100):
            yax[i]=vdb.M(z[i],M0[j])/M0[j]
            yax2[i]=M(z[i],M0[j])/M0[j]
        plt.plot(np.log10(1+z),np.log10(yax),'--',color=color)
        plt.plot(np.log10(1+z),np.log10(yax2),color=color,label='5e'+str(j+11)+'$M_\odot$')

    data1=np.loadtxt('5e11.dat')
    data1=np.transpose(data1)
    data2=np.loadtxt('5e12.dat')
    data2=np.transpose(data2)
    data3=np.loadtxt('5e13.dat')
    data3=np.transpose(data3)
    data4=np.loadtxt('5e14.dat')
    data4=np.transpose(data4)

    plt.plot(np.log10(1.+data1[1]),data1[3],':',color='k')
    plt.plot(np.log10(1.+data2[1]),data2[3],':',color='r')
    plt.plot(np.log10(1.+data3[1]),data3[3],':',color='b')
    plt.plot(np.log10(1.+data4[1]),data4[3],':',color='g')

    plt.legend()    
    plt.show()
    
    # z=np.linspace(0,6,100)
    # M0=1.e14*Msun
    # yax=np.zeros(100)
    # for i in range(100):
    #     yax[i]=M(z[i],M0)/Msun
    #     print(z[i],'%3g'%yax[i])

if __name__ == "__main__":
    main()

