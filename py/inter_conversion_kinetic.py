# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library
from matplotlib.widgets import Slider
from scipy.integrate import odeint

__package__ = "1st Order Kinetic reaction"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"


def conversion (y, t, kf, kb):

  product, reagents = y
  dp =  kf * reagents - kb * product
  dr = -kf * reagents + kb * product

  return (dp, dr)

def update (val):

  kf = skf.val
  kb = skb.val
  res = odeint(conversion, y0, time, args=(kf, kb))

  plast, rlast = res[-1]
  peq = get_peq(plast, rlast, kf, kb)
  req = get_req(plast, rlast, kf, kb)

  if not np.isclose(peq, plast, rtol=1e-3) or not np.isclose(req, rlast, rtol=1e-3):
    raise ValueError('Unstable Dynamic')

  l1.set_ydata(res[:, 0])
  l2.set_ydata(res[:, 1])
  l3.set_ydata(np.sum(res, axis=1))
  fig.canvas.draw_idle()


get_peq = lambda p, r, kf, kb : (kf * (p + r)) / (kf + kb)
get_req = lambda p, r, kf, kb : (kb * (p + r)) / (kf + kb)


if __name__ == '__main__':


  '''
  A first order kinetic is given by

          R -><- P

        dR/dt = -kf*R + kb*P
        dP/dt = kf*R - kb*P

        R(t) + P(t) = Cost
  '''

  # initial conditions
  y0 = (0, 1.)
  kf, kb = (.3, .6)

  time = np.linspace(0, 10, 1000)
  res = odeint(conversion, y0, time, args=(kf, kb))
  total = np.sum(res, axis=1)

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
  fig.subplots_adjust(bottom=0.25)

  ax.set_xlabel('Time', fontsize=14)
  ax.set_ylabel('', fontsize=14)
  ax.set_title('1st Order Kinetic Reaction (Inter-conversion)', fontsize=14)

  l1, l2,  = ax.plot(time, res)
  l1.set_label('Products')
  l2.set_label('Reagents')
  l3, = ax.plot(time, total, 'g-', label='total', alpha=.5, linestyle='dashed')

  ax.margins(x=0)

  axcolor = 'lightgoldenrodyellow'
  axkf = plt.axes([.25, .02, .65, .03], facecolor=axcolor)
  axkb = plt.axes([.25, .07, .65, .03], facecolor=axcolor)

  skf    = Slider(axkf, r'kf', 1e-2, 1., valinit=kf, valstep=.1)
  skb    = Slider(axkb, r'kb', 1e-2, 1., valinit=kb, valstep=.1)

  skf.on_changed(update)
  skb.on_changed(update)
  ax.legend(loc='best', fontsize=14)

  plt.show()
