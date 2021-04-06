from matplotlib import pyplot as plt 
import numpy as np 
data=np.loadtxt('analysis.dat')
data=np.transpose(data)
unitrho=1.67262171e-27
unitl=3.0856775807e21
unitv=1.e7
unitmdot=unitrho*unitl**2*unitv
msun=2.e33
yr=3.145e7
unitmdot=unitmdot/msun*yr
unittime=unitl/unitv
unittime=unittime/yr/1.e9
# plt.semilogy(data[0]*unittime,-data[3]*unitmdot,label='Cold gas')
plt.semilogy(data[0]*unittime,-data[2]*unitmdot,label='All gas')
plt.ylim(1e-1,1e4)
plt.title("Infalling gas")
plt.xlabel("Time (Gyr)")
plt.ylabel("Mdot ($M_\odot/yr$)")
plt.legend()
plt.show()
