import numpy as np 
from matplotlib import pyplot as plt 

data=np.loadtxt('cool.txt')
data=np.transpose(data)
plt.semilogx(data[0],data[1])
plt.show()
# print(data)