#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library

__package__ = "Newton Method"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

f = lambda x : x**2 - 2 # function
d = lambda x : 2*x # derivative
newton = lambda x : x - f(x) / d(x) # newton formula

if __name__ == '__main__':

  n = 10000
  toll = 1e-4
  x0 = 30
  x = np.linspace(start=0, stop=10, num=n)
  y = f(x)

  step = [x0]
  tst = 1
  while tst >= toll:
  	step.append( newton(x0) ) # array of steps (redundant)
  	x0 = step[-1] # update root
  	tst = abs(step[-1] - step[-2]) # check convergence

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.plot(x, y, 'b-', label="f(x)") # plot function
  ax.plot(step, f(np.asarray(step)), 'r--', label="newton step") # plot steps
  ax.plot(step[-1], f(step[-1]), 'go') # plot root
  ax.set_xlabel("x", fontsize=14)
  ax.set_ylabel("f(x)", fontsize=14)
  ax.set_title("Newton method", fontsize=14)

  plt.show()
