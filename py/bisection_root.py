#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library

__package__ = "Bisection method"
__author__  = "Nico Curti"
__email__   = "nico.curit2@unibo.it"

f = lambda x : x**2 - 2
bisection = lambda x0, x1 : (x0, x1) if f(x0)*f(x1) < 0 else (x1, x0)

if __name__ == '__main__':

  n = 10000
  toll = 1e-4
  x0 = 5
  x1 = 1
  x = np.linspace(start = 0, stop = 10, num = n)
  y = f(x)
  step = [x0]
  tst = 1
  while tst >= toll:
    c = (x0 - x1) * .5 # middle point
    x0, x1 = bisection(c, x0) # bisection update
    step.append(c) # array of steps
    tst = abs(step[-1] - step[-2]) # check convergence

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.plot(x, y, 'b-', label="f(x)") # plot function
  ax.plot(step, f(np.asarray(step)), 'r--', label="bisection step") # plot steps
  ax.plot(step[-1], f(step[-1]), 'go') # plot root
  ax.set_xlabel("x", fontsize=14)
  ax.set_ylabel("f(x)", fontsize=14)
  ax.set_title("Bisection method", fontsize=14)
