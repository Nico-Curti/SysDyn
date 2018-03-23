# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:30:58 2018

@author: NICO
"""

import numpy as np
import matplotlib.pylab as plt
import scipy

dt = .01
it = 10000
sigma = 16.
b = 4.
r = 45.92
x, y, z = np.empty(shape=(it)), np.empty(shape=(it)), np.empty(shape=(it))
x[0], y[0], z[0] = 10, 1, 1

Vx = lambda x, y, sigma : sigma*(y-x)
Vy = lambda x, y, z, r : -x*z + r*x - y
Vz = lambda x, y, z, b : -b*z + x*y

for i in range(0, it-1):
	x[i+1] = x[i] + dt*Vx(x[i], y[i], sigma)
	y[i+1] = y[i] + dt*Vy(x[i], y[i], z[i], r)
	z[i+1] = z[i] + dt*Vz(x[i], y[i], z[i], b)

# False Nearest Neighbours

RT = 15
AT = 2
sigmay = np.std(x)

maxEmbDim = 10
delay= 1
rEEM = it - (maxEmbDim*delay - delay)
EEM = np.concatenate([y[delay*maxEmbDim-delay:].reshape(1, len(y[delay*maxEmbDim-delay:])), 
                      [y[delay*maxEmbDim-(i+1)*delay:-i*delay] for i in range(1, maxEmbDim)]
                      ], axis=0).T
Ind1 = np.empty(maxEmbDim)
Ind2 = np.empty(maxEmbDim)
embedm = 0 # only for plot
for k in range(1, maxEmbDim+1):
    D = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(EEM[:, :k], "euclidean"))        
    np.fill_diagonal(a=D, val=np.inf)
    l = np.argmin(D[:rEEM - maxEmbDim - k, :], axis=1)
    fnn1 = np.asarray([abs(y[i + maxEmbDim + k - 1]-y[li + maxEmbDim + k - 1])/D[i, li] for i,li in enumerate(l) if D[i, li] > 0 and li + maxEmbDim + k - 1 < it])
    fnn2 = np.asarray([abs(y[i + maxEmbDim + k - 1]-y[li + maxEmbDim + k - 1])/sigmay for i,li in enumerate(l) if D[i, li] > 0 and li + maxEmbDim + k - 1 < it])
    Ind1[k-1] = len(np.where(np.asarray(fnn1) > RT)[0])
    Ind2[k-1] = len(np.where(np.asarray(fnn2) > AT)[0])
    
    if embedm == 0: # only for plot
	    if Ind1[k-1] / len(fnn1) < .1 and Ind2[-1] / len(fnn1) < .1 and Ind1[k-1] != 0:
	        embedm = k
	        #break # uncomment for true algorithm
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
ax.plot(np.arange(0, maxEmbDim), Ind1)
ax.set_xlabel('Embedding dimension', fontsize=14)
ax.set_ylabel('% FNN', fontsize=14)
ax.set_title('Optimal Embedding Dimension with FNN', fontsize=16)
ax.plot(embedm, Ind1[embedm], 'r.') 
plt.text(embedm, Ind1[embedm] + 100, "EmbDim = $%d$"%(embedm))
