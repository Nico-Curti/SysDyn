#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library
import pandas as pd            # Dataframe
import seaborn as sns          # advance graphics
from scipy.integrate import odeint # python integrator

__package__ = "Toggle Switch"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

def toggle(state, t, αx, αy, βx, βy):
  x, y = state
  dx = αx/(1+y**βx) - x
  dy = αy/(1+x**βy) - y
  return (dx, dy)

def integrate_toggle(state_0, time, args):
  states = odeint(toggle, state_0, time, args)
  results = pd.DataFrame(states, columns=['x', 'y'])
  results['time'] = time
  return results

if __name__ == '__main__':

  initial_state = (0.5, 0.5)
  αx = 1.1
  αy = 0.9
  βx = βy = 1.1
  params = (αx, αy, βx, βy)
  time = np.linspace(0, 10, 101)
  stati = odeint(toggle, initial_state, time, args=params)
  results = pd.DataFrame(stati, columns=['x', 'y'])
  results['time'] = time

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
  ax.plot('time', 'x', data=results)
  ax.plot('time', 'y', data=results)
  ax.legend()
  sns.despine(fig, trim=True, offset=10)

#%%

  fig, ax = plt.subplots(figsize=(12, 8))

  αx = αy = 2
  βx = βy = 3
  params = (αx, αy, βx, βy)
  for x in np.linspace(0, 2, 100):
    y = 2-x
    state_0 = [x, y]
    time = np.linspace(0, 10, 101)
    results = integrate_toggle(state_0, time, params)
    ax.plot('time', 'x', data=results,
            label=str(x),
            color='b',
            alpha=0.5)
  fig.legend(loc='best', fontsize=14)

#%%

  fig, axis = plt.subplots(figsize=(12, 8))

  αx = αy = 2

  for beta in np.linspace(1, 3, 101):
    for x in np.linspace(0, 2, 10):
      y = 2-x
      state_0 = [x, y]
      time = np.linspace(0, 100, 2)
      params = (αx, αy, beta, beta)
      results = integrate_toggle(state_0, time, params)
      last_x = results['x'].values[-1]
      axis.scatter(beta, last_x,
              color='b',
              alpha=0.5)

  #%%

  res = {}
  for alpha in [1, 2, 3]:
    for beta in np.linspace(1, 3, 50):
      for x in np.linspace(0, 2, 10):
        state_0 = [x, 2-x]
        time = [0, 100]
        params = (alpha, alpha, beta, beta)
        results = odeint(toggle, state_0, time, params)
        last_x = results[-1, 0]
        res[alpha, beta, x] = last_x

  data = pd.Series(res).reset_index()
  data.columns = ['alpha', 'beta', 'start_x', 'last_x']

  #%%

  fg = sns.FacetGrid(data, row='alpha')
  fg.map(plt.scatter, 'beta', 'last_x')

  #%%

  from itertools import product

  alpha_x_seq = [1, 2, 3]
  alpha_y_seq = [1, 2, 3]
  beta_seq = np.linspace(1, 4, 800)
  x_seq = [0, 2]

  to_iter = product(alpha_x_seq, alpha_y_seq, beta_seq, x_seq)
  res = dict()
  for (alpha_x, alpha_y, beta, x) in to_iter:
    state_0 = [x, 2-x]
    time = [0, 100]
    params = (alpha_x, alpha_y, beta, beta)
    results = odeint(toggle, state_0, time, params)
    last_x = results[-1, 0]
    res[alpha_x, alpha_y, beta, x] = last_x

  data = pd.Series(res).reset_index()
  data.columns = ['alpha_x', 'alpha_y', 'beta',
                  'start_x', 'last_x']

  # %%

  fg = sns.FacetGrid(data, row='alpha_x', col='alpha_y')
  fg.map(plt.scatter, 'beta', 'last_x')


