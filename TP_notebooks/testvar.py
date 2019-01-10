from burgers import *
from gausscov import *
from simvar import *
from obsopt import *
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
sigmab = 0.01              # background state error std
sigmao = 0.001             # Observation error std
Lb = 0.05                  # Correlation length for B matrix

# Assimilation Parameters

precond = True              # preconditioning by square root of B 
iobstsub = 5                # Frequency of temporal subsampling of observations, [1:nt], 1=every time step
iobsxsub = 4                # Frequency of spatial subsampling of observations, [1:nx], 1=every space step

# Observation operator and error covariance matrix

H = Obsopt(nx,iobsxsub,nt,iobstsub)
R = sigmao*sigmao*np.eye(nx,nx)

# Initialization of true field uo
uo=np.sin(2*math.pi*xx);
true=H.gen_obs(M,uo,sigmao)

# Initialization of background
ub=np.cos(2*math.pi*xx)
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

var=Variational(ubkg[0],nt,B,M,H,R,precond)

J=0.
alpha=0.001

uopt= np.random.normal(0.,sigmab,uo.size)

Jini  = var.cost(uopt)
grini = var.grad(uopt)
norm  = grini.dot(grini)
for iii in range(1,21):
    uctl = uopt + alpha * grini
    J = var.cost(uctl)
    print alpha, (J-Jini)/(alpha*norm)
    alpha /= 10.

