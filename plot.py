import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time as time
import plot

print plot.metropolis_hastings.__doc__

t = time.time()

X = np.array([range(25),]*25)
Y = np.array([range(25),]*25).transpose()
bucket = plot.metropolis_hastings()

print time.time() - t

user_input = raw_input("choose a layer (0 to 24) ")
while user_input.isdigit():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(X, Y, bucket[int(user_input),:,:])
    plt.show()

    user_input = raw_input("choose a layer (0 to 24)")

