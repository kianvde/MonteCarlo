import h2_mh
import numpy as np
import helper
import matplotlib.pyplot as plt

# print the h2_mh documentation
h2_mh_file = open('h2_mh.f90', 'r')
for line in h2_mh_file:
    if (line[0] != '!'):
        break
    print line,

s = 0.5
beta = 0.7
a = helper.find_a(s)

h2_mh.init_rnd()
r = 4.0*np.random.rand(6)-2.0
energy, descent, r = h2_mh.metropolis_hastings(r, [s, a, beta], 0.2, 200000, False)
b_avg = np.array(float("Inf"))
E_avg = np.array(float("Inf"))
crit = float("Inf")
i = 0
while crit > (0.01 / 100.0):
    b_temp = np.array([0.0])
    E_temp = np.array([0.0])
    for i2 in range(10):
        energy, descent, r = h2_mh.metropolis_hastings(r, [s, a, beta], 0.2, 5000, True)
        beta -= 0.1*descent
        b_temp = np.hstack((b_temp,beta))
        E_temp = np.hstack((E_temp,energy))
    avg_temp = np.sum(b_temp[1:])/10.0
    Eavg_temp = np.sum(E_temp[1:])/10.0
    b_avg = np.hstack((b_avg,avg_temp))
    E_avg = np.hstack((E_avg,Eavg_temp))
    crit = abs(b_avg[i+1] - b_avg[i])
    i += 1

    print beta

# energy, descent, r = h2_mh.metropolis_hastings(r, [s, a, beta], 0.2, 1000000, True)
# print "energy: "
# print energy
# print "beta: "
# print beta

b_calc = np.array(float("Inf"))
E_calc = np.array(float("Inf"))
for k in range(20):
    b_temp = np.array([0.0])
    E_temp = np.array([0.0])
    for k2 in range(10):
        energy, descent, r = h2_mh.metropolis_hastings(r, [s, a, beta], 0.2, 5000, True)
        beta -= 0.1*descent
        b_temp = np.hstack((b_temp,beta))
        E_temp = np.hstack((E_temp,energy))
    avg_temp = np.sum(b_temp[1:])/10.0
    Eavg_temp = np.sum(E_temp[1:])/10.0
    b_calc = np.hstack((b_calc,avg_temp))
    E_calc = np.hstack((E_calc,Eavg_temp))
    print b_calc[k+1]
b_hope = np.average(b_calc[1:])
E_hope = np.average(E_calc[1:])
print b_hope

b_plot = np.hstack((b_avg[1:],b_calc[1:]))
avg_plot = np.ones((1,np.size(b_plot)))*b_hope
E_plot = np.hstack((E_avg[1:],E_calc[1:]))
Eavg_plot = np.ones((1,np.size(E_plot)))*E_hope

plt.figure()
plt.plot(range(np.size(E_plot)),np.transpose(E_plot))
plt.plot(range(np.size(E_plot)),np.transpose(Eavg_plot))
plt.plot(range(np.size(b_plot)),np.transpose(b_plot))
plt.plot(range(np.size(b_plot)),np.transpose(avg_plot))
plt.show()


    # print "energy: "
    # print energy
    # print "beta: "
    # print beta

