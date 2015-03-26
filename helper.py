import numpy as np
import math
import h2_mh


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


def find_steady_state(s,a,beta,r):
    energy, descent, r = h2_mh.metropolis_hastings(r, [s, a, beta], 1.0, 200000, False)
    return r


def minimize_beta(s,a,beta,eps,r):
    avg_old = float("Inf")
    i = 0
    crit = float("Inf")
    while crit > (eps / 10.0 / 100.0):
        b_temp = np.array([0.0])
        E_temp = np.array([0.0])
        for i in range(10):
            energy, descent, r = h2_mh.metropolis_hastings(r, [s, a, beta], 1.0, 50000, True)
            beta -= 0.1*descent
            b_temp = np.hstack((b_temp,beta))
            E_temp = np.hstack((E_temp,energy))
        crit = abs(np.sum(b_temp[1:])/10.0 - avg_old)
        avg_old = np.sum(b_temp[1:])/10.0
    return beta, r









