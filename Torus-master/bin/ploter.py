import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.gca(projection='3d')

r,pr,phi,pphi,z,pz,t= np.loadtxt('R.txt', delimiter=',', unpack=True,skiprows=1)
x=r*np.cos(phi)
y=r*np.sin(phi)
ax.plot(x, y, z, label='orbit')
ax.legend()
plt.show()
