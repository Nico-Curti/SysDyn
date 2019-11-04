# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library
from matplotlib.widgets import Slider
from scipy.integrate import odeint

__package__ = "1st Order Kinetic reaction"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"


def conversion (y, t, alpha):

  product, reagents = y
  dp =  alpha * reagents
  dr = -alpha * reagents

  return (dp, dr)

def update (val):

  alpha = salpha.val
  y0 = (sp0.val, sr0.val)
  res = odeint(conversion, y0, time, args=(alpha,))
  l1.set_ydata(res[:, 0])
  l2.set_ydata(res[:, 1])
  l3.set_ydata(np.sum(res, axis=1))
  fig.canvas.draw_idle()


if __name__ == '__main__':


  '''
  A first order kinetic is given by

          dR/dt = - k*R
          dP/dt =   k*R

  like the previous one equation. In fact we have already proved that the
  solution is an exponential decay.
  '''

  # initial conditions
  y0 = (.5, .5)
  alpha = .5

  time = np.linspace(0, 10, 1000)
  res = odeint(conversion, y0, time, args=(alpha,))
  total = np.sum(res, axis=1)

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
  fig.subplots_adjust(bottom=0.25)

  ax.set_xlabel('Time', fontsize=14)
  ax.set_ylabel('', fontsize=14)
  ax.set_title('1st Order Kinetic Reaction', fontsize=14)

  l1, l2,  = ax.plot(time, res)
  l1.set_label('Products')
  l2.set_label('Reagents')
  l3, = ax.plot(time, total, 'g-', label='total', alpha=.5, linestyle='dashed')

  ax.margins(x=0)

  axcolor = 'lightgoldenrodyellow'
  axalpha = plt.axes([.25, .02, .65, .03], facecolor=axcolor)
  axp0    = plt.axes([.25, .07, .65, .03], facecolor=axcolor)
  axr0    = plt.axes([.25, .12, .65, .03], facecolor=axcolor)

  salpha = Slider(axalpha, r'$\alpha$', 0., 1.5, valinit=alpha, valstep=.1)
  sp0    = Slider(axp0,    r'p0', 0., 1., valinit=y0[0], valstep=.1)
  sr0    = Slider(axr0,    r'r0', 0., 1., valinit=y0[1], valstep=.1)

  salpha.on_changed(update)
  sp0.on_changed(update)
  sr0.on_changed(update)
  ax.legend(loc='best', fontsize=14)

  plt.show()
