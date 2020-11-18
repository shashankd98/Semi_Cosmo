import numpy as np 
from matplotlib import pyplot as plt 
from scipy.optimize import newton

sec_year=3.154e7
h=0.72
H0=100*1e3/3.086e22*h
omega_m=0.27
omega_l=1.-omega_m
Msun=2.e33
M0=1e11*Msun
zf=4.
nu=1.211+1.858*np.log10(1+zf)+0.308*omega_l**2-0.032*np.log10(M0/(1e11/h*Msun))

def t(z):
    time=1/H0*2/3/np.sqrt(omega_l)*np.log((np.sqrt(omega_l/(1+z)**3)+np.sqrt(omega_l/(1+z)**3+omega_m))/np.sqrt(omega_m))
    return time

def t_inv(z,t0):
    return t(0)-t(z)-t0

def M(z):
    rhs=-0.301*(np.log10(1+z)/np.log10(1+zf))**nu
    return M0*10**(rhs)

def main():
    t0=t(0)
    z=np.linspace(0,10,100)
    Mz=M(z)
    # print(M(zf)/M0)

    z4=10**((np.log10(0.04)/-0.301)**(1/nu)*np.log10(1+zf))-1
    # print((t0-t(z4))/3.154e7/1e9)
    t_test=1e10*sec_year
    z_test=newton(t_inv,1.5,args=[t_test])
    # print(z_test)
    # tz=t(z)
    plt.plot(z,(t(0)-t(z))/1e9/3.154e7)
    plt.show()
    # print(t(1.5)/sec_year/1e9)
    # plt.plot(np.log10(1.+z),np.log10(Mz/M0))
    # plt.show()
if __name__ == "__main__":
    main()