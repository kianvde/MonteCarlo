import h2_mh
import numpy as np

# print the h2_mh documentation
h2_mh_file = open('h2_mh.f90', 'r')
for line in h2_mh_file:
    if (line[0] != '!'):
        break
    print line,

beta = 1.0
for i in range(5):
    r = 4.0*np.random.rand(6)-2.0
    h2_mh.init_rnd()
    energy, descent, r = h2_mh.metropolis_hastings(r, [2.0, 0.90182909, beta], 1.0, 5000, False)
    energy, descent, r = h2_mh.metropolis_hastings(r, [2.0, 0.90182909, beta], 1.0, 10000000, True)
    beta -= 10*descent

    print "energy: "
    print energy
    print "beta: "
    print beta

