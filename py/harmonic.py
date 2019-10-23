#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np             # numerical library
import matplotlib.pylab as plt # plot library

__package__ = "Harmonic oscillator"
__author__  = "Nico Curti"
__email__   = "nico.curit2@unibo.it"

Vx = lambda q, p : p
Vy = lambda q, p : -q

# Integrators
euler        = lambda q, p, dt, k : ( q + dt*Vx(q, p), p + dt*k*Vy(q, p) )
simplettic   = lambda q, p, dt, k : ( q + dt*Vx(q, p + dt*k*Vy(q, p)), p + dt*k*Vy(q, p))

if __name__ == '__main__':
  k  = .5
  it = 1000
  dt = .1
  q_eu, p_eu = np.empty(it), np.empty(it)
  q_eu[0], p_eu[0] = .5, 0 # initial conditions
  q_si, p_si = np.empty(it), np.empty(it)
  q_si[0], p_si[0] = .5, 0 # initial conditions

  for t in range(0, it-1):
    q_eu[t+1], p_eu[t+1] = euler(q_eu[t], p_eu[t], dt, k)
    q_si[t+1], p_si[t+1] = simplettic(q_si[t], p_si[t], dt, k)

  time = np.arange(0, it*dt, dt)
  fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
  ax1.plot(time, q_eu, 'b-', label="q euler")
  ax1.plot(time, p_eu, 'r-', label="p euler")
  ax2.plot(q_eu, p_eu, 'b-', label="phase space")

  ax3.plot(time, q_si, 'b-', label="q simplettic")
  ax3.plot(time, p_si, 'r-', label="p simplettic")
  ax4.plot(q_si, p_si, 'b-', label="phase space")

  plt.show()
