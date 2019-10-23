#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np                            # numerical library
import matplotlib.pylab as plt                # plot library
from sklearn.metrics import mutual_info_score # mutual information algorithm
from mpl_toolkits.mplot3d import Axes3D       # 3D plot

__package__ = "Mutual information"
__author__  = "Nico Curti"
__email__   = "(nico.curit2@unibo.it)"

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
  x[0], y[0], z[0] = (10, 1, 1) # initial conditions

  # euler integration
  for i in range(0, it-1):
    x[i+1] = x[i] + dt*Vx(x[i], y[i], sigma)
    y[i+1] = y[i] + dt*Vy(x[i], y[i], z[i], r)
    z[i+1] = z[i] + dt*Vz(x[i], y[i], z[i], b)


  # Optimal time delay with Mutual information
  mutual_info = lambda x, y, bins : mutual_info_score(None,
                            None,
                            contingency=np.histogram2d(x, y, bins)[0])
  tau_max = 100
  mi = np.empty(shape=(tau_max-1,), dtype=float)
  for tau in range(1, tau_max):
    # mutual information between delay signal and original
    mi[tau-1] = mutual_info(x[tau:], x[:-tau], 100)
  optimal_delay = np.where(np.diff(mi) > 0)[0][0] # find first minimum

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  ax.plot(np.arange(1, tau_max)*dt, mi, 'bo-')
  ax.set_xlabel("delay", fontsize=14)
  ax.set_ylabel("mutual information", fontsize=14)

  ax.scatter(dt*np.arange(1, tau_max)[optimal_delay],
             mi[optimal_delay],
             marker = 'o',
             s=400,
             facecolors="none",
             edgecolor="r")

  # Phase space reconstruction
  fig = plt.figure(figsize=(8,8))
  ax = Axes3D(fig)
  ax.plot(x[:-optimal_delay*2], x[optimal_delay:-optimal_delay], x[:-2*optimal_delay])

  plt.show()
