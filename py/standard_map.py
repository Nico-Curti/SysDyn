#!/usr/bin/python

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library

__package__ = "Pendulum Standard map vectorized"
__author__  = "Nico Curti (nico.curit2@unibo.it)"

if __name__ == '__main__':
  it    = 100
  size  = 5000
  k     = 1.
  Q     = np.empty(shape=(size, it))
  P     = np.empty(shape=(size, it))
  Q[0]  = np.random.rand(it)
  P[0]  = np.random.rand(it)
  pi_2  = 2. * np.pi
  inv_2pi = 1. / pi_2

  for i in range(size-1):
    P[i+1] = np.mod(P[i] - (k * inv_2pi * np.sin(pi_2 * Q[i])), 1)
    Q[i+1] = np.mod(Q[i] + P[i+1] + .5, 1)

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.plot(Q, P, 'k,')
  ax.set_xlim(0, 1)
  ax.set_ylim(0, 1)
