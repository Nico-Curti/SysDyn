#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np                       # numerical library
import matplotlib.pylab as plt           # plot library
import matplotlib.animation as animation # animation plot

__package__ = "Diffusion 2D model"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

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

if __name__ == '__main__':

  A     = 4.5
  B     = 6.75
  Du    = 2.
  Dv    = 16.
  dt    = .005
  dx    = 1
  Tmax  = 50000
  dim   = 120
  delta = .3

  ut = A   + delta*np.random.rand(dim, dim)
  vt = B/A + delta*np.random.rand(dim, dim)

  ims = np.empty(Tmax + 1, dtype=np.object)

  # FFMpegWriter = animation.writers['ffmpeg']
  # writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)

  fig = plt.figure(figsize=(8,8))
  sizes = ut.shape
  fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward=False)
  ax = plt.Axes(fig, [0., 0., 1., 1.])
  ax.set_axis_off()
  fig.add_axes(ax)

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
  # movie.save('diffusion.mp4', writer=writer)
  plt.show()
