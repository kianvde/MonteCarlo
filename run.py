import h2_mh
import numpy as np
import helper
import matplotlib.pyplot as plt

# # print the h2_mh documentation
# h2_mh_file = open('h2_mh.f90', 'r')
# for line in h2_mh_file:
#     if (line[0] != '!'):
#         break
#     print line,

s_vector = np.linspace(1,2,20)
energy = np.zeros(np.size(s_vector))
for i, s in enumerate(s_vector):
    # problem parameters
    beta = 0.7
    sigma = 0.2

    a = helper.find_a(s)
    parameters = [s, a, beta]

    # initialize fortran random value generator and r
    h2_mh.init_rnd()
    r = 4.0*np.random.rand(6)-2.0

    # get mh updater to steady state
    sigma, parameters[2], r = helper.find_beta(parameters, sigma, r)

    # calculate the energy
    energy[i], descent, acceptance, r = h2_mh.metropolis_hastings(r, parameters, sigma, 1000000)

    print "hallo"


plt.figure()
plt.plot(s_vector, energy)
plt.show()