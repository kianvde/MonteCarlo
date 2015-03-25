# calculate the energy for the hydrogen molecule using
# the metropolis hastings algorithm

import time as time
import h2_mh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

t = time.time()

# parameter 2 (a): is a trial wave function parameter
# parameter 3 (beta): is a trial wave function parameter
# parameter 1 (s): is the seperation between the protons
#                   following from a*(1+exp(-s/a))=1
energy, plot_bin = h2_mh.metropolis_hastings(2.0, 0.90182909, 0.001)
print energy

print time.time() - t

X = np.array([np.linspace(-2,2,50),]*50)
Y = np.array([np.linspace(-2,2,50),]*50).transpose()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, plot_bin)
plt.show()