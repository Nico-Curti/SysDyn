#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np                      # numerical library
import matplotlib.pylab as plt          # plot library
from mpl_toolkits.mplot3d import Axes3D # 3D plot

__package__ = "Lorentz attractor"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

# lorentz formula
Vx = lambda x, y, sigma : sigma*(y-x)
Vy = lambda x, y, z, r : -x*z + r*x - y
Vz = lambda x, y, z, b : -b*z + x*y

if __name__ == '__main__':

  dt      = .01
  it      = 10000
  sigma   = 16.
  b       = 4.
  r       = 45.92
  x, y, z = np.empty(shape=(it), dtype=float), np.empty(shape=(it), dtype=float), np.empty(shape=(it), dtype=float)
  x[0], y[0], z[0] = 10, 1, 1 # initial conditions

  # euler integration
  for i in range(0, it-1):
    x[i+1] = x[i] + dt*Vx(x[i], y[i], sigma)
    y[i+1] = y[i] + dt*Vy(x[i], y[i], z[i], r)
    z[i+1] = z[i] + dt*Vz(x[i], y[i], z[i], b)

  fig = plt.figure(figsize=(8,8))
  ax  = Axes3D(fig)
  ax.plot(x, y, z)

  plt.show()
