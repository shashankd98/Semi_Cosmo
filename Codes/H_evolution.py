import numpy as np 

h=0.7
z=6
Ol=0.6911
Om=0.3089

hz=h*(Om*(1+z)**3+Ol)**(0.5)
print(hz)
