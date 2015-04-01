import numpy as np
import math
import h2_mh
import time

# find a from s using Newton-Rhapson
def find_a(s):
    
    a = s
    while True:
        f = 1.0/(1.0 + math.exp(-s/a)) - a

        f_prime = -s*math.exp(-s/a)/(a*math.exp(-s/a))**2 - 1.0
        a -= f / f_prime

        if abs(f / f_prime) < 10**-10:
            return a

# find beta using mh algorithm
# returns sigma, beta, r
def find_beta(parameters, sigma, r):
    # update parameters
    gamma = 0.1     # damped gradient descent factor
    delta = 0.05    # acceptance rate update factor

    for i in range(500):
        energy, descent, acceptance, r = h2_mh.metropolis_hastings(r, parameters, sigma, 5000)
        sigma += delta*(acceptance - 0.5)
        parameters[2] -= gamma*descent

    return parameters[2], sigma, r

# function for the Morse potential
def morse(x, a, D, r0, c):
    return D*(1 - np.exp(-a*(x-r0)))**2 - c










