from gausscov import *
from analyseKF import *
from burgers import *
from obsopt import *

import numpy as np
import math



# Space-time domain
nx = 40                     # number of grid points
dx = 1./nx                  # space step
xx = np.array(range(nx))*dx # grid points abscissa
dt = 0.5*dx                 # time step
nt = 20                     # number of time steps
ns = 0                      # numerical scheme

M=Burgers(nx,dx,dt,ns)

# Error staristics
sigmab = 0.001              # background state error std
sigmao = 0.001             # Observation error std
Lb = 0.05                  # Correlation length for B matrix

# Assimilation Parameters

iobstsub = 1                # Frequency of temporal subsampling of observations, [1:nt], 1=every time step
iobsxsub = 10                # Frequency of spatial subsampling of observations, [1:nx], 1=every space step

# Observation operator and error covariance matrix

H = Obsopt(nx,iobsxsub,nt,iobstsub)
R = sigmao*sigmao*np.eye(H.nobs,H.nobs)

# Initialization of true field uo
uo=np.sin(2*math.pi*xx);
yo=np.dot(H.mat,uo) + np.random.normal(0.,sigmao,H.nobs)
                                       
# Initialization of background
ub=np.cos(2*math.pi*xx)

# Initialization of Pf matrix and its sqare root
    
B = gausscov(nx,sigmab,Lb,2)
Pf = B.mat
Sf = B.sqr

ua,Sa = analyseKF(ub,Sf,H.mat,yo,R)

Pa = np.dot(Sa , Sa.T)

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

print np.dot(H.mat,xx)

f, axarr = plt.subplots(2, 2)

axarr[0, 0].plot(xx,uo,'k-')
axarr[0, 0].plot(xx,ub,'b-')
axarr[0, 0].plot(xx,ua,'r-',linewidth=3)
axarr[0, 0].plot(np.dot(H.mat,xx),yo,'kd')
axarr[0, 0].legend(['True','Background','Analysis','Observations'])
axarr[0, 0].set_title('BLUE analysis')

axarr[0, 1].set_title('BLUE increment')
axarr[0, 1].plot(xx,ua-ub,'m-',linewidth=3)
axarr[0, 1].plot(np.dot(H.mat,xx),np.zeros(H.nobs),'kd')
axarr[0, 1].legend(['Increment','Observations'])

cmap = plt.get_cmap('PiYG')
levels = MaxNLocator(nbins=15).tick_values(-Pf.max(), Pf.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
axarr[1, 0].pcolormesh(xx, xx, Pf,cmap=cmap, norm=norm)
axarr[1, 0].set_title('Pf')

axarr[1, 1].pcolormesh(xx, xx, Pa,cmap=cmap, norm=norm)
axarr[1, 1].set_title('Pa')



plt.show()


