#!/usr/bin/python

import numpy as np # numerical library
import matplotlib.pylab as plt # plot library
from scipy.integrate import odeint # python integrator

__package__ = "Logistic Map"
__author__  = "Nico Curti (nico.curit2@unibo.it)"

logistic = lambda x, t, mu : mu * (np.mod(x, 1) - np.mod(x**2, 1)) # logistic formula

if __name__ == '__main__':

  it   = 60
  dt   = .1
  x0   = np.mod(.1, 1)
  mu   = 1.3
  t    = np.arange(0, dt*it, dt)
  x    = np.empty(shape=(it))
  x[0] = x0

  for i in range(0, it-1):
    x[i+1] = np.mod( np.mod(x[i], 1) + dt*logistic(x[i], t[i], mu) , 1)

  xt = odeint(logistic, x0, t, args=(mu, ))

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.scatter(t, x, marker='o', facecolors="none", edgecolor='k', s=20, label="euler")
  ax.plot(t, xt, 'r-', lw=2, alpha=.5, label="scipy")
  ax.set_ylim(0, 1)
  ax.legend(loc='best', fontsize=14)
