# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:30:58 2018

@author: NICO
"""

import numpy as np # numerical library
import matplotlib.pylab as plt # plot library
import matplotlib.animation as animation # animation plot

#%% Bernoulli2D
bernoulli2D = lambda x : np.mod(2*x, 1) # bernoulli formula

# initial condition
N = 100 # number of points
x = np.linspace(0, 1, N) # x steps
y = np.linspace(0, 1, N) # y steps
meanx, stdx = np.mean(x), np.std(x)
meany, stdy = np.mean(y), np.std(y)
x, y = np.meshgrid(x, y)
G = np.exp( - ( .5*(x-meanx)**2 / stdx**2 + .5*(y-meany)**2 / stdy**2 ) )


fig = plt.figure(figsize=(8,8))
time = 100 # number of iterations
ims = np.empty(time - 1, dtype=np.object)
ev0 = G
for i in range(1, time):
    ev = bernoulli2D(ev0)
    ims[i-1] = [plt.imshow( ev, animated=True )]
    ev0 = ev
movie = animation.ArtistAnimation(fig,
        								  ims,
                                   interval=50, 
                                   blit=True,
                                   repeat_delay=100
                                  )