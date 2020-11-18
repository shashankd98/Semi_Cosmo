import numpy as np 
from matplotlib import pyplot as plt 

Msun=2.e33

# data=np.loadtxt("res.txt")
# data=np.transpose(data)
# fig = plt.figure(figsize=(20,10))
# ax1=fig.add_subplot(221)
# ax2=fig.add_subplot(222)
# ax3=fig.add_subplot(223)
# ax4=fig.add_subplot(224)

# ax1.loglog(data[0],data[1])
# ax2.loglog(data[0],data[2])
# ax3.loglog(data[0],data[3])
# ax4.semilogx(data[0],data[4])
# ax3.loglog(data[0],data[5])
# plt.show()

data=np.loadtxt("m.txt")
data=np.transpose(data)
z=data[0]
M=data[1]
plt.plot(np.log10(1+z),np.log10(M/M[0]))
plt.show()