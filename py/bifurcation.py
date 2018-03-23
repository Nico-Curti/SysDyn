# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:30:58 2018

@author: NICO
"""

import numpy as np # numerical library
import matplotlib.pylab as plt # plot library

# logistic map
logistic = lambda x, mu : 4*mu*x*(1-x)

# Bifurcation logistic map x_{n+1} = \mu x_n (1 - x_n)
it = 100000
mu = np.linspace(.5, 1, it) # coefficient values
x = np.empty(it)
x[0] = .5

for t in range(0, it-1):
	x[t+1] = logistic(x[t], mu[i])

critical_y = x[np.where(np.diff(x) > 0)[0][0]]
critical_x = mu[critical_y]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
ax.scatter(mu, x, marker=',', color='k')
ax.vlines(critical_x, 0, 1, color='r', linestyle="dashed")
ax.set_xlabel("$\mu$", fontsize=14)
ax.set_ylabel("$x_k$", fontsize=14)
ax.text(critical_y, critical_x, "$\mu = %.3f$"critical[0])
ax.set_xlim(.5, 1.01)