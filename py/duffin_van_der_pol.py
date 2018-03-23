# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:30:58 2018

@author: NICO
"""

import numpy as np # numerical library
import matplotlib.pylab as plt # plot library
from mpl_toolkits.mplot3d import Axes3D # 3D plot

# Duffin-Van Der Pol oscillator
dv_der_pol = lambda x, y, z, mu, f, o, dt: (dt*y, dt*(mu*1-x*x)*y - x*x*x + f*np.cos(z), z + dt*o)

mu    = .2
f     = 1.
omega = .94
dt    = 1e-3
it    = 10000

x, y, z = np.empty(it), np.empty(it), np.empty(it)
x[0], y[0], z[0] = 1., 0., 1. # initial conditions
time = np.arange(0, dt*it, it) # time

for t in range(0, it-1):
	x[t+1], y[t+1], z[t+1] = dv_der_pol(x[t], y[t], z[t], mu, f, o, dt)

fig = plt.figure(figsize=(8,8))
ax = Axes3D(fig)
ax.plot(x, y, np.sin(z))
