# -*- coding: utf-8 -*-
#!/usr/bin/env python

import numpy as np                   # numerical library
import matplotlib.pylab as plt       # plot library
from scipy.integrate import odeint   # ODE integrator
from scipy.optimize import curve_fit # fit function

__package__ = "Fit example"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

def derivate (y, t, alpha=1):
  return -alpha * y

def decay (x, alpha):

  y = odeint(derivate, y0, t, args=(alpha,))
  return y[:, 0]


if __name__ == '__main__':

  # initial conditions
  y0 = 10
  alpha = 0.5

  # time intervals
  t = np.linspace(start=0, stop=5, num=100)

  # integrate the ODE (aka the model)
  y = odeint(derivate, y0, t, args=(alpha,))

  # generate synthetic data adding some noise

  y_obs = y + np.random.normal(loc=0., scale=.5, size=y.shape)

  params, covariances = curve_fit(f=decay,
                                  xdata=t,
                                  ydata=y_obs[:, 0])

  y_est = decay(t, params)

  print('Fit Error: {:.3f}'.format(np.sqrt(covariances[0][0])))

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
  ax.scatter(t, y_obs, marker='o', alpha=.5, label='data points')
  ax.plot(t, y_est, color='blue', label='fit')
  ax.plot(t, y, color='red', linestyle='dashed', label='theory')
  ax.set_xlabel(r'$x$', fontsize=14)
  ax.set_ylabel(r'$f(x)$', fontsize=14)
  ax.legend(loc='upper right', fontsize=14)

  plt.show()
