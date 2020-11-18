import numpy as np 
from matplotlib import pyplot as plt  

x=np.linspace(0,9,1000)
beta=4
gamma=8
f=(1+(x)**beta)**(-gamma/beta)
t=np.polyfit(x,f,3)
fit=np.poly1d(t)
tx=fit(x)

plt.plot(x,f)
# plt.plot(x,tx)
plt.show()