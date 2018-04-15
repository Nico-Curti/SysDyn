# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 16:48:28 2018

@author: NICO
"""

import numpy as np # numerical library
import matplotlib.pylab as plt # plot library
import matplotlib.animation as animation # animation plot

g = lambda u, v, A, B : B*u - v*u*u
f = lambda u, v, A, B : A - (B+1)*u + v*u*u

def laplacian(X, dx):
    n, m = X.shape
    lap = np.empty((n, m))
    lap[1:-1, 1:-1] = (X[1:-1, 2:] + X[1:-1, :-2] + X[2:, 1:-1] + X[:-2, 1:-1] - 4*X[1:-1, 1:-1]) / dx / dx
    # Laplacian with boundary conditions
    lap[1:-1, 0] = (X[1:-1, 1] + X[1:-1, -1] + X[:-2, 0] + X[2:, 0] - 4*X[1:-1, 0]) / dx / dx 
    lap[1:-1,-1] = (X[1:-1, 0] + X[1:-1, -2] + X[:-2,-1] + X[2:,-1] - 4*X[1:-1,-1]) / dx / dx 
    lap[0, 1:-1] = (X[0, 2:  ] + X[0, :-2  ] + X[-1,1:-1]+ X[1,1:-1] - 4*X[0, 1:-1]) / dx / dx 
    lap[-1,1:-1] = (X[-1,2:  ] + X[-1,:-2  ] + X[-2,1:-1]+ X[0,1:-1] - 4*X[-1,1:-1]) / dx / dx 
    lap[0, 0]  = (X[0, 1] + X[1, 0] + X[0,-1] + X[-1, 0] - 4*X[0, 0]) / dx / dx
    lap[-1,-1] = (X[-1,-2]+ X[-2,-1]+ X[-1,0] + X[0, -1] - 4*X[-1,-1])/ dx / dx
    lap[0, -1] = (X[0,-2] + X[1,-1] + X[0, 0] + X[-1, 0] - 4*X[0,-1]) / dx / dx
    lap[-1, 0] = (X[-1,1] + X[0, 0] + X[-1,-1]+ X[-2, 0] - 4*X[-1,0]) / dx / dx
    return lap

A = 4.5
B = 6.75
Du = 2.
Dv = 16.
dt = .005
dx = 1
Tmax = 50000
dim = 120
delta = .3

ut = A   + delta*np.random.rand(dim, dim)
vt = B/A + delta*np.random.rand(dim, dim)

ims = np.empty(Tmax + 1, dtype=np.object)
fig = plt.figure(figsize=(8,8))
ims[0] = [plt.imshow(ut, animated=True, cmap="jet")]
for t in range(Tmax):
    u, v = ut, vt
    ut = dt * (Du * laplacian(u, dx) + f(u, v, A, B)) + u
    vt = dt * (Dv * laplacian(v, dx) + g(u, v, A, B)) + v
    ims[t + 1] = [plt.imshow( ut, animated=True, cmap="jet" )]
movie = animation.ArtistAnimation(fig,
                                  ims,
                                  interval=50, 
                                  blit=True,
                                  repeat_delay=100
                                  )