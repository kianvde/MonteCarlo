# calculate the energy for the hydrogen molecule using
# the metropolis hastings algorithm

import time as time
import h2_mh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# parameter 1 (s): is the separation between the protons
# parameter 2 (a): is a trial wave function parameter
#                  following from a*(1+exp(-s/a))=1
# parameter 3 (beta): is a trial wave function parameter
beta = 1.0

for i in range(5):
    energy, descent = h2_mh.metropolis_hastings(2.0, 0.90182909, beta)
    beta -= 10*descent
    print "energy: "
    print energy
    print "beta: "
    print beta

# plt.figure()
# plt.plot(energy)
# plt.show()

# plot
# X = np.array([np.linspace(-2,2,50),]*50)
# Y = np.array([np.linspace(-2,2,50),]*50).transpose()
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_wireframe(X, Y, plot_bin)
# plt.show()