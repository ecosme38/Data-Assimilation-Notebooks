from burgers import *
import math
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Space-time domain
nx = 40                    # number of grid points
dx = 1./nx                  # space step
xx = np.array(range(nx))*dx          # grid points abscissa
dt = 0.5*dx                # time step
nt = 20                    # number of time steps
ns=0                       # numerical scheme

model=Burgers(nx,dx,dt,ns)
# Initialization of field uu
uu=np.sin(2.*math.pi*xx)
umat=list()
umat.append(uu)                     # Storage for future plot

# Reference trajectory

for it in range(nt):
    uu=model.step(uu)
    umat.append(uu)

# plot

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'r-', animated=True)

def init():
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, 1)
    return ln,

def update(frame):
    xdata=xx #xdata.append(xx)
    ydata=frame #ydata.append(frame)
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=umat,
                    init_func=init, blit=True)
plt.show()
