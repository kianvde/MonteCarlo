import h2_mh
import numpy as np
from scipy.optimize import curve_fit
import helper
import matplotlib.pyplot as plt

# problem parameters
beta = 0.7
sigma = 0.2
s = np.linspace(1,2,50)

energy = np.zeros(np.size(s))
for i,s_i in enumerate(s):

    # find a from a*(1-e(-s/a)) = 1 using the Newton Rhapson method
    a = helper.find_a(s_i)
    parameters = [s_i, a, beta]

    # initialize fortran random value generator and r
    h2_mh.init_rnd()
    r = 4.0*np.random.rand(6) - 2.0

    # find the correct value for beta
    parameters[2], sigma, r = helper.find_beta(parameters, sigma, r)

    # calculate the energy
    energy[i], descent, acceptance, r = h2_mh.metropolis_hastings(r, parameters, sigma, 1000000)

    print i+1

# fit the morse potential
popt, pcov = curve_fit(helper.morse, s, energy)
x = np.linspace(1,2,1000)
var = np.sqrt(np.diag(pcov))

print "minimum energy: " + "{:.5f}".format(popt[3]) + " +/-" + "{:.5f}".format(var[3])
print "D: " + "{:.5f}".format(popt[1]) + " +/-" + "{:.5f}".format(var[1])
print "r0: " + "{:.5f}".format(popt[2]) + " +/-" + "{:.5f}".format(var[2])

# plot the fit and data
plt.figure()
plt.plot(x, helper.morse(x, popt[0], popt[1], popt[2], popt[3]))
plt.plot(s, energy, 'x')
plt.show()


