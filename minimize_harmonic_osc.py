# this script minimizes the wavefunction sqrt(pi/a)*exp(-a*x^2) for
# the harmonic potential V=.5*x^2
#
# with hbar = omega = m = 1 theory predicts 1/2 for a

import harm_osc as ho
import numpy as np
import matplotlib.pyplot as plt

a = np.linspace(.3, .7, 30)
energy = np.zeros(np.size(a))
for i in range(np.size(a)):
    energy[i] = ho.metropolis_hastings(a[i])

plt.figure()
plt.plot(a,energy)
plt.show()



