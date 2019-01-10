from burgers import *
from gausscov import *
from simvar import *
from obsopt import *
from plots import *

import numpy as np
import scipy.optimize as opt
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
sigmab = 0.02              # background state error std
sigmao = 0.001             # Observation error std
Lb = 0.05                  # Correlation length for B matrix

# Assimilation Parameters

precond = True             # preconditioning by square root of B (1=yes)
iobstsub = 5                # Frequency of temporal subsampling of observations, [1:nt], 1=every time step
iobsxsub = 8                # Frequency of spatial subsampling of observations, [1:nx], 1=every space step

# Observation operator and error covariance matrix

H = Obsopt(nx,iobsxsub,nt,iobstsub)
R = sigmao*sigmao*np.eye(H.nobs,H.nobs)


# Initialization of true field uo
uo=np.sin(2*math.pi*xx);
true=H.gen_obs(M,uo,sigmao)

# Initialization of background
ub=np.cos(2.*math.pi*xx)
ubkg=[ub]
for it in range(nt):
    ub=M.step(ub)
    ubkg.append(ub)

# Initialization of B matrix and its inverse

if precond:
    indic=2
else:
    indic=1
    
B=gausscov(nx,sigmab,Lb,indic)

# Actual minimisation

var=Variational(ubkg[0],nt,B,M,H,R,precond)

#print uo-uopt
#err=opt.check_grad(var.cost,var.grad,uo,epsilon=1.e-15)
#print err

#approx=opt.approx_fprime(uo,var.cost,1.e-15)
#print approx
#print var.grad(uo)

uopt=np.zeros(nx)

res = opt.minimize(var.cost,uopt,
                   method='L-BFGS-B',
                   jac=var.grad,
                   options={'disp': True, 'gtol': 1e-05, 'maxiter': 10000, 'iprint':100})

print (res)


if precond:
    ua=ubkg[0] + B.sqr.dot(res['x'])
else:
    ua=ubkg[0] + res['x']

uana=[ua]
for it in range(nt):
    ua=M.step(ua)
    uana.append(ua)

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

f, axarr = plt.subplots(2, 2)

axarr[0, 0].plot(xx,true[0],'k-')
axarr[0, 0].plot(xx,ubkg[0],'b-')
axarr[0, 0].plot(xx,uana[0],'r-',linewidth=3)
axarr[0, 0].plot(np.dot(H.mat,xx),H.yo[0],'kd')
axarr[0, 0].legend(['True','Background','Analysis','Observations'])
axarr[0, 0].set_title('States at the begining of experiments')

axarr[0, 1].plot(xx,true[nt],'k-')
axarr[0, 1].plot(xx,ubkg[nt],'b-')
axarr[0, 1].plot(xx,uana[nt],'r-',linewidth=3)
axarr[0, 1].plot(np.dot(H.mat,xx),H.yo[nt],'kd')
axarr[0, 1].legend(['True','Background','Analysis','Observations'])
axarr[0, 1].set_title('States at the end of experiments')

axarr[1, 0].set_title('Errors at the begining of experiments')
axarr[1, 0].plot(xx,(ubkg[0]-true[0])**2,'b-',linewidth=3)
axarr[1, 0].plot(xx,(uana[0]-true[0])**2,'r-',linewidth=3)
axarr[1, 0].legend(['Squared background error','Squared analysis error'])
    
axarr[1, 1].set_title('Errors at the end of experiments')
axarr[1, 1].plot(xx,(ubkg[nt]-true[nt])**2,'b-',linewidth=3)
axarr[1, 1].plot(xx,(uana[nt]-true[nt])**2,'r-',linewidth=3)
axarr[1, 1].legend(['Squared background error','Squared analysis error'])
    
plt.show()

# fig, ax = plt.subplots()
# ax.pcolormesh(xx, xx,np.eye(nx,nx)-res['hess_inv'].todense())
# plt.show()

anim(xx,nt,[true,ubkg,uana],legends=['True','Background','Analysis'])
plt.show

