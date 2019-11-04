#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import scipy
import pandas as pd
import seaborn as sns

__package__ = "Grassberger-Procaccia algorithm"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

Vx = lambda x, y, sigma : sigma*(y-x)
Vy = lambda x, y, z, r : -x*z + r*x - y
Vz = lambda x, y, z, b : -b*z + x*y


if __name__ == '__main__':
  dt      = .01
  it      = 10000
  sigma   = 16.
  b       = 4.
  r       = 45.92
  x, y, z = np.empty(shape=(it)), np.empty(shape=(it)), np.empty(shape=(it))
  x[0], y[0], z[0] = 10, 1, 1

  for i in range(0, it-1):
    x[i+1] = x[i] + dt*Vx(x[i], y[i], sigma)
    y[i+1] = y[i] + dt*Vy(x[i], y[i], z[i], r)
    z[i+1] = z[i] + dt*Vz(x[i], y[i], z[i], b)

  # Grassberger Procaccia

  omit_pts = 3
  ED = scipy.spatial.distance.pdist(np.concatenate([x.reshape(len(x), 1), x.reshape(len(x),1)]), "euclidean")
  min_eps = ED[ED!=0].min()
  max_eps = 2**np.ceil(np.log(ED.max())/np.log(2))
  eps_vec_size = int(np.floor( (np.log(max_eps/min_eps) / np.log(2)))) + 1
  eps_vec = max_eps*2.**(-np.arange(0, eps_vec_size))

  Npairs = it*(it-1)*.5
  C_eps = np.asarray([len(ED[(ED < eps) & (ED > 0)])/Npairs for eps in eps_vec])

  k1 = omit_pts + 1
  k2 = eps_vec_size - omit_pts
  xp = np.log(eps_vec[np.arange(k1, k2+1)]) / np.log(2)
  yp = np.log(C_eps[np.arange(k1, k2+1)]) / np.log(2)

  data = pd.DataFrame(data=[xp, yp], index=["xp", "yp"]).T
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
  plot = sns.regplot(x="xp",
                     y="yp",
                     data=data,
                     ax=ax).legend()

  plt.show()
