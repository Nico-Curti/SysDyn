# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:30:58 2018

@author: NICO
"""

import numpy as np # numerical library
import matplotlib.pylab as plt # plot library

#%% Pendulum Standard map vectorized

Q = np.empty(shape=(int(it/10), it))
P = np.empty(shape=(int(it/10), it))
Q[0] = np.random.rand(it)
P[0] = np.random.rand(it)
for i in range(int(it/10)-1):
    P[i+1] = np.mod(P[i] + dt*(k**2*(-np.sin(Q[i]))), 1)
    Q[i+1] = np.mod(Q[i] + dt*P[i+1], 1)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
ax.plot(Q, P, 'k,')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)