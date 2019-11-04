# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library
from matplotlib.widgets import Slider
from scipy.integrate import odeint

__package__ = "Michaelis Menten"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"


def conversion (y, t, kf, kb, Kirr):

  S, E, ES, P = y
  d_S  = -kf * S * E + kb * ES
  d_E  = -kf * S * E + kb * ES + Kirr * ES
  d_ES =  kf * S * E - kb * ES - Kirr * ES
  d_P  = Kirr * ES
  return (d_S, d_E, d_ES, d_P)

def update (val):

  kf = skf.val
  kb = skb.val
  Kirr = ski.val
  res = odeint(conversion, y0, time, args=(kf, kb, Kirr,))
  l1.set_ydata(res[:, 0])
  l2.set_ydata(res[:, 1])
  l3.set_ydata(res[:, 2])
  l4.set_ydata(res[:, 3])
  fig.canvas.draw_idle()


if __name__ == '__main__':

  # initial conditions
  y0 = (10, 1, 0, 0) # (S0, E0, ES0, P0)
  kf, kb, Kirr = (1., 1e-2, 1.)
  params = (kf, kb, Kirr)

  time = np.geomspace(1e-2, 20, 101)
  res = odeint(conversion, y0, time, args=params)

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
  fig.subplots_adjust(bottom=0.25)

  ax.set_xlabel('Time', fontsize=14)
  ax.set_ylabel('', fontsize=14)
  ax.set_title('Michaelis Menten', fontsize=14)

  l1, l2, l3, l4  = ax.plot(time, res)
  l1.set_label('$dS$')
  l2.set_label('$dE$')
  l3.set_label('$dES$')
  l4.set_label('$dP$')

  ax.margins(x=0)

  axcolor = 'lightgoldenrodyellow'
  axkf = plt.axes([.25, .02, .65, .03], facecolor=axcolor)
  axkb = plt.axes([.25, .07, .65, .03], facecolor=axcolor)
  axki = plt.axes([.25, .12, .65, .03], facecolor=axcolor)

  skf = Slider(axkf, r'$k_f$', 0., 10., valinit=kf, valstep=.1)
  skb = Slider(axkb, r'$k_b$', 0., 10., valinit=kb, valstep=.1)
  ski = Slider(axki, r'Kirr',  0., 10., valinit=Kirr, valstep=.1)

  skf.on_changed(update)
  skb.on_changed(update)
  ski.on_changed(update)
  ax.legend(loc='best', fontsize=14)

  plt.show()
