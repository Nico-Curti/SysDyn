#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library

__package__ = "Bifurcation logistic map x_{n+1} = \\mu x_n (1 - x_n)"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

logistic = lambda x, mu : 4*mu*x*(1-x)


if __name__ == '__main__':

  it   = 100000
  mu   = np.linspace(.5, 1, it) # coefficient values
  x    = np.empty(it)
  x[0] = .5

  for t in range(0, it-1):
    x[t+1] = logistic(x[t], mu[t])

  critical_y = x[np.argmax(np.diff(x) > 0)]
  critical_x = mu[int(critical_y)]

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.scatter(mu, x, marker=',', color='k', s=1)
  ax.vlines(critical_x, 0, 1, color='r', linestyle="dashed")
  ax.set_xlabel("$\\mu$", fontsize=14)
  ax.set_ylabel("$x_k$", fontsize=14)
  ax.text(critical_y+0.1, critical_x, "$\\mu = %.3f$"%critical_x)
  ax.set_xlim(.5, 1.01)

  plt.show()
