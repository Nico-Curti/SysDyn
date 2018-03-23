# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:30:58 2018

@author: NICO
"""

import numpy as np # numerical library

#%% Integration

f = lambda x : np.sin(x)
rettangoli = lambda dx, N : sum( [dx * f(i * dx) for i in range(N)] )
trapezi = lambda dx, N : sum( [dx * (f(i*dx) + f((i+1)*dx)) *.5 for i in range(N)] )
simpson = lambda dx, N : sum( [dx/6*(f(i*dx) + 4*f((i*dx + (i+1)*dx)*.5) + f((i+1)*dx)) for i in range(N)] )

truth = 1. - np.cos(1.) # sin primitive in [0, 1]
dx = 1.
for n in range(1, 100):
	N = n*n
	print( N, abs(rettangoli(dx / N, N) - truth), abs(trapezi(dx / N, N) - truth), abs(simpson(dx / N, N) - truth) )

