import numpy as np
import math
import h2_mh

# find a from s using Newton Rhapson
def find_a(s):
    a = s
    crit = float("Inf")
    while crit > 10**-10:
        f = 1.0/(1.0 + math.exp(-s/a)) - a

        fprime_pre = 1.0/(a*math.exp(-s/a))**2
        fprime = -s*math.exp(-s/a)*fprime_pre - 1.0

        a -= f / fprime
        crit = abs(f / fprime)

    return a

# find beta using mh algorithm
# returns sigma, beta, r
def find_beta(parameters, sigma, r):

    # update parameters
    gamma = 0.05    # damped gradient descent factor
    delta = 0.05    # acceptance rate update factor

    for i in range(1000):
        energy, descent, acceptance, r = h2_mh.metropolis_hastings(r, parameters, sigma, 5000)
        sigma += delta*(acceptance - 0.5)
        parameters[2] -= gamma*descent
    return sigma, parameters[2], r









