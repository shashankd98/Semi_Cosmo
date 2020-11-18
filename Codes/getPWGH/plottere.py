import numpy as np 
from matplotlib import pyplot as plt 

data1=np.loadtxt('5e11.dat')
data1=np.transpose(data1)
data2=np.loadtxt('5e12.dat')
data2=np.transpose(data2)
data3=np.loadtxt('5e13.dat')
data3=np.transpose(data3)
data4=np.loadtxt('5e14.dat')
data4=np.transpose(data4)

plt.plot(np.log10(1.+data1[1]),data1[3])
plt.plot(np.log10(1.+data2[1]),data2[3])
plt.plot(np.log10(1.+data3[1]),data3[3])
plt.plot(np.log10(1.+data4[1]),data4[3])
plt.show()