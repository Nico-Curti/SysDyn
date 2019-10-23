#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library

__package__ = "Secant Method"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

f = lambda x: x**2 - 3 # function
secant = lambda x0, x1 : x0 - (x1 - x0) / (f(x1) - f(x0)) * f(x0) # secant formula

if __name__ == '__main__':

  n    = 10000
  toll = 1e-4
  x0   = 2
  x1   = 0
  x    = np.linspace(start=0, stop=10, num=n)
  y    = f(x)
  step = [x0]
  tst  = 1

  while tst >= toll:
    step.append( secant(x0, x1) ) # array of steps (redundant)
    x0 = x1 # update x0
    x1 = step[-1] # update x1
    tst = abs(step[-1] - step[-2]) # check convergence

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.plot(x, y, 'b-', label="f(x)") # plot function
  ax.plot(step, f(np.asarray(step)), 'r--', label="secant step") # plot step
  ax.plot(step[-1], f(step[-1]), 'go') # plot root
  ax.set_xlabel("x", fontsize=14)
  ax.set_ylabel("f(x)", fontsize=14)
  ax.set_title("Secant method", fontsize=14)

  plt.show()
