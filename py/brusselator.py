# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library
from matplotlib.widgets import Slider

__package__ = "Brusselator"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"


Vx = lambda x, y, A, B : A + x*x*y - B*x - x
Vy = lambda x, y, A, B : B*x - x*x*y


def RK4 (x, y, dt, A, B):

  kx1 = dt * Vx(x, y, A, B)
  ky1 = dt * Vy(x, y, A, B)

  kx2 = dt * Vx(x * .5 * kx1, y + .5 * ky1, A, B)
  ky2 = dt * Vy(x * .5 * kx1, y + .5 * ky1, A, B)

  kx3 = dt * Vx(x + .5 * kx2, y + .5 * ky2, A, B)
  ky3 = dt * Vy(x + .5 * kx2, y + .5 * ky2, A, B)

  kx4 = dt * Vx(x + kx3, y + ky3, A, B)
  ky4 = dt * Vy(x + kx3, y + ky3, A, B)

  x_new = x + 1. / 6 * (kx1 + 2. * kx2 + 2. * kx3 + kx4)
  y_new = y + 1. / 6 * (ky1 + 2. * ky2 + 2. * ky3 + ky4)

  return (x_new, y_new)


def update (val):

  A = sA.val
  B = sB.val

  for i in range(0, iterations):
    x[i + 1], y[i + 1] = RK4(x[i], y[i], dt, A, B)

  l1.set_ydata(x[:-1])
  l2.set_ydata(y[:-1])
  l3.set_xdata(x)
  l3.set_ydata(y)

  ax1.set_ylim(min(min(x), min(y)), max(max(x), max(y)))

  ax2.set_xlim(min(min(x), min(x)), max(max(x), max(x)))
  ax2.set_ylim(min(min(y), min(y)), max(max(y), max(y)))

  fig.canvas.draw_idle()



if __name__ == '__main__':

  # initial conditions
  A, B = (.5, 2.)

  dt = 1e-2
  iterations = 10000

  time = np.linspace(0, dt * iterations, iterations)
  x = np.empty(shape=(iterations + 1), dtype=float)
  y = np.empty(shape=(iterations + 1), dtype=float)

  x[0] = 1.6
  y[0] = 2.8

  for i in range(0, iterations):
    x[i + 1], y[i + 1] = RK4(x[i], y[i], dt, A, B)
  # 2.46 s ± 145 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

  fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
  fig.subplots_adjust(bottom=0.25)

  ax1.set_xlabel('Time', fontsize=14)
  ax1.set_ylabel('', fontsize=14)
  ax1.set_title('Brusselator', fontsize=14)

  l1, = ax1.plot(time, x[:-1], color='r', label='x')
  l2, = ax1.plot(time, y[:-1], color='b', label='y')

  ax2.set_xlabel('X', fontsize=14)
  ax2.set_ylabel('Y', fontsize=14)
  ax2.set_title('Phase Space', fontsize=14)

  l3, = ax2.plot(x, y, color='b')

  ax1.margins(x=0)
  ax2.margins(x=0)

  axcolor = 'lightgoldenrodyellow'
  axA = plt.axes([.25, .07, .65, .03], facecolor=axcolor)
  axB = plt.axes([.25, .12, .65, .03], facecolor=axcolor)

  sA = Slider(axA, r'A', 0., 8., valinit=A, valstep=.1)
  sB = Slider(axB, r'B', 0., 8., valinit=B, valstep=.1)

  sA.on_changed(update)
  sB.on_changed(update)
  ax1.legend(loc='best', fontsize=14)

  plt.show()
