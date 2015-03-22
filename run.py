# calculate the energy for the hydrogen molecule using
# the metropolis hastings algorithm

import time as time
import h2_mh

t = time.time()

# parameter 1 (a): is a trial wave function parameter
# parameter 2 (beta): is a trial wave function parameter
# parameter 3 (s): is the seperation between the protons
#                   following from a*(1+exp(-s/a))=1
energy = h2_mh.metropolis_hastings(2, 0.90182909, 1.0)
print energy

print time.time() - t

