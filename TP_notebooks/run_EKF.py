from gausscov import *
from analyseKF import *
from burgers import *
from obsopt import *
from plots import *

import numpy as np
import math

# Space-time domain
nx = 40                     # number of grid points
dx = 1./nx                  # space step
xx = np.array(range(nx))*dx # grid points abscissa
dt = 0.5*dx                 # time step
nt = 40                     # number of time steps
ns = 0                      # numerical scheme

M=Burgers(nx,dx,dt,ns)

# Error staristics
sigmab = 0.01              # background state error std
sigmao = 0.01             # Observation error std
Lb = 0.05                  # Correlation length for B matrix

# Assimilation Parameters

iobstsub = 5                # Frequency of temporal subsampling of observations, [1:nt], 1=every time step
iobsxsub = 8                # Frequency of spatial subsampling of observations, [1:nx], 1=every space step

# Observation operator and error covariance matrix

H = Obsopt(nx,iobsxsub,nt,iobstsub)
R = sigmao*sigmao*np.eye(H.nobs,H.nobs)

# Initialization of true field uo and true trajectory
uo=np.sin(2*math.pi*xx);
true=H.gen_obs(M,uo,sigmao)

# Initialization of background
ub=np.cos(2*math.pi*xx)
ubkg=[ub]
for it in range(nt):
    ub=M.step(ub)
    ubkg.append(ub)

# Initialization of Pf matrix and its sqare root
    
B = gausscov(nx,sigmab,Lb,2)
P = B.mat
S = B.sqr
uu=ubkg[0]
uana=[]
ufor=[uu]
Pfmat=[np.real(np.diag(P))]
Pamat=[]

#------------  KALMAN FILTER   ----------------------
# -------------------------------------------------------
# --- TIME LOOP ---
for it in range(nt):

    # -- ANALYSIS ---------------------------------
    
    if H.isobserved(it):
  
        up=uu
        Sp=S
        uu,S=analyseKF(up,Sp,H.mat,H.yo[it],R)
        
    uana.append(uu)
    Pamat.append(np.real(np.diag(np.dot(S,S.T))))

  
    # -- FORECAST ------------------------------
    # Mean state
    up=uu
    uu=M.step(up)
    # Error modes (square root of cov. matrix)
    for imem in range(nx):
        uerrp = up + S[:,imem]
        uerr = M.step(uerrp)
        S[:,imem]=uerr-uu
  
    ufor.append(uu)
    Pfmat.append(np.real(np.diag(np.dot(S,S.T))))
    
# --- END OF TIME LOOP ---
# Last analysis, if obs exists after the last time step:
if H.isobserved(nt):
    up=uu
    Sp=S
    uu,S=analyseKF(up,Sp,H.mat,H.yo[nt],R)

uana.append(uu)
Pamat.append(np.real(np.diag(np.dot(S,S.T))))

P=np.dot(S,S.T) # For P diagnostics if desired



import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

f, axarr = plt.subplots(2, 2)

axarr[0, 0].plot(xx,true[nt],'k-')
axarr[0, 0].plot(xx,ubkg[nt],'b-')
axarr[0, 0].plot(xx,ufor[nt],'g*')
axarr[0, 0].plot(xx,uana[nt],'r-',linewidth=3)
axarr[0, 0].plot(np.dot(H.mat,xx),H.yo[nt],'kd')
axarr[0, 0].legend(['True','Background','Forecast','Analysis','Observations'])
axarr[0, 0].set_title('States at end of experiment')

axarr[0, 1].set_title('End of experiment')
axarr[0, 1].plot(xx,(uana[nt]-true[nt])**2,'m-',linewidth=3)
axarr[0, 1].plot(xx,Pamat[nt]/(sigmab*sigmab),'r-')
axarr[0, 1].legend(['Squared analysis error','P$^a$ variance(rescaled by $\sigma_b^2$)'])

rmse=[]
tme=[]
for i in range(nt+1):
    errf=ufor[i]-true[i]
    erra=uana[i]-true[i]
    rmse.append(np.mean(errf*errf))
    tme.append(i)
    if np.array_equal(errf,erra):
        rmse.append(np.mean(errf*errf))
        tme.append(i)
        

axarr[1, 0].set_title('RMS error vs time')
axarr[1, 0].plot(tme,rmse,'g-',linewidth=3)

cmap = plt.get_cmap('PiYG')
levels = MaxNLocator(nbins=15).tick_values(-np.real(P.max()), np.real(P.max()))
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
axarr[1, 1].pcolormesh(xx, xx, P,cmap=cmap, norm=norm)
axarr[1, 1].set_title('P$^a$')

plt.show()

# Animations

animation(xx,nt,[true,ubkg,uana,ufor],legends=['True','Background','Analysis','Forecast'])

for i in range(nt+1):
    true[i]=uana[i]-true[i]
    ubkg[i]=uana[i]-ubkg[i]
    ufor[i]=(uana[i]-ufor[i])
    Pamat[i]=Pamat[i]/(sigmab*sigmab)

animation(xx,nt,[true,ubkg],legends=['Analysis-reference','Analysis-background'])
animation(xx,nt,[ufor,Pamat],legends=['Analysis-forecast','Analysis variance(rescaled by $\sigma_b^2$)'])

