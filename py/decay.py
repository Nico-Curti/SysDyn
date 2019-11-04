# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library
from matplotlib.widgets import Slider

__package__ = "0 Order Kinetic reaction"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"


def update (val):

  alpha = salpha.val
  y0 = sy0.val
  y = integrate(x, y0, alpha)
  l.set_ydata(y)
  ax.set_ylim(0, y0)
  fig.canvas.draw_idle()


def integrate (x, y0, alpha):

  dx = x[1] - x[0]
  y = np.empty(shape=(len(x) + 1, ), dtype=float)
  y[0] = y0

  for i, _ in enumerate(x):
    y[i + 1] = y[i] * (1. - alpha * dx)

  return y[:-1]

if __name__ == '__main__':


  '''
  A 0 order kinetic equation is given by

              R -> P

  You should already know its solution but
  know you can also evaluate this formula
  using a numerical integration

          dR/dt = - k

          (R_{i+1} - R_{i})/dt = -k
          R_{i+1} = -k * dt + R_{i}

  and we can integrate it.

  In this example we see a slight more complex
  version given by the decay equation and we try
  to integrate it using the Euler method
  '''

  # initial conditions
  y0 = 10
  alpha = .5
  iterations = 1000
  dt = 1e-2

  x = np.arange(0, iterations*dt, dt)
  y = integrate(x, y0, alpha)

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
  fig.subplots_adjust(bottom=0.25)

  ax.set_ylabel(r'$y$', fontsize=14)
  ax.set_xlabel(r'$x$', fontsize=14)
  ax.set_title('Decay', fontsize=14)

  l, = ax.plot(x, y, label=r'$\hat{y} = - \alpha \cdot y$')

  ax.margins(x=0)

  axcolor = 'lightgoldenrodyellow'
  axalpha = plt.axes([.25, .02, .65, .03], facecolor=axcolor)
  axy0    = plt.axes([.25, .07, .65, .03], facecolor=axcolor)

  salpha = Slider(axalpha, r'$\alpha$', 0., 1.5, valinit=alpha, valstep=.1)
  sy0    = Slider(axy0,    r'y0', 0., 20, valinit=y0, valstep=1.)

  salpha.on_changed(update)
  sy0.on_changed(update)
  ax.legend(loc='best', fontsize=14)

  plt.show()
