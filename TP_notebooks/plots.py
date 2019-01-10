import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class animator:
    ''' Setup init and update methods for matplotlib.animation.FuncAnimation.
    not very nicely done, but I can't make it work properly with variable-size trajectories
    '''
    def __init__(self,xx,ax=None,xmin=0,xmax=1,ymin=-1,ymax=1,trajectories=None,legends=None,colors=['k-','b-','r-','g-']):
        self.xx=xx
        self.xmax=xmax
        self.xmin=xmin
        self.ymin=ymin
        self.ymax=ymax
        self.ax=ax
        if not isinstance(trajectories[0], (list, tuple)) :
            trajectories=[trajectories]
        self.trajectories=trajectories
        self.ncurve=len(trajectories)
        if self.ncurve>4:
            raise ValueError('Cannot animate more than 4 curves at once')
        self.legends=legends
        self.colors=colors
        self.xdata, ydata = [], []
        
    def init(self):
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)
        if self.legends!=None:
            self.ax.legend(self.legends)
            
        if self.ncurve==1:
            self.ln1, = self.ax.plot([], [], self.colors[0], animated=True)
            return self.ln1,
        elif self.ncurve==2:
            self.ln1, = self.ax.plot([], [], self.colors[0], animated=True)
            self.ln2, = self.ax.plot([], [], self.colors[1], animated=True)
            return self.ln1,self.ln2,
        elif self.ncurve==3:
            self.ln1, = self.ax.plot([], [], self.colors[0], animated=True)
            self.ln2, = self.ax.plot([], [], self.colors[1], animated=True)
            self.ln3, = self.ax.plot([], [], self.colors[2], animated=True)
            return self.ln1,self.ln2,self.ln3,
        elif self.ncurve==4:
            self.ln1, = self.ax.plot([], [], self.colors[0], animated=True)
            self.ln2, = self.ax.plot([], [], self.colors[1], animated=True)
            self.ln3, = self.ax.plot([], [], self.colors[2], animated=True)
            self.ln4, = self.ax.plot([], [], self.colors[3], animated=True)
            return self.ln1,self.ln2,self.ln3,self.ln4,
        

    def update(self,i):
        self.xdata=self.xx 
        self.ydata=self.trajectories[0][i]
        self.ln1.set_data(self.xdata,self.ydata)
        if self.ncurve>1:
            self.ydata=self.trajectories[1][i]
            self.ln2.set_data(self.xdata,self.ydata)
        if self.ncurve>2:
            self.ydata=self.trajectories[2][i]
            self.ln3.set_data(self.xdata,self.ydata)
        if self.ncurve>3:
            self.ydata=self.trajectories[3][i]
            self.ln4.set_data(self.xdata,self.ydata)
            
        if self.ncurve==1:
            return self.ln1
        elif self.ncurve==2:
            return self.ln1,self.ln2,
        elif self.ncurve==3:
            return self.ln1,self.ln2,self.ln3,
        elif self.ncurve==4:
            return self.ln1,self.ln2,self.ln3,self.ln4,


def anim(xx, nt, trajectories,**kwargs):

    if not isinstance(trajectories[0], (list, tuple)) :
        trajectories=[trajectories]

    if len(trajectories)==1 :
        animation_1(xx, nt, trajectories, **kwargs)
    elif len(trajectories)==2 :
        animation_2(xx, nt, trajectories, **kwargs)
    elif len(trajectories)==3 :
        animation_3(xx, nt, trajectories, **kwargs)
    elif len(trajectories)==4 :
        animation_4(xx, nt, trajectories, **kwargs)

        
def animation_1(xx, nt, trajectories,legends=None,colors=['k-']):

    if isnotebook():
        plt.rcParams["animation.html"] = "jshtml"

    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln1, = ax.plot([], [], colors[0], animated=True)
    
    def init():
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)
        if legends!=None:
            ax.legend(legends)
        return ln1,

    def update(frame):
        xdata=xx 
        ydata=frame
        ln1.set_data(xdata,ydata)
        return ln1,

    ani = FuncAnimation(fig, update, frames=self.trajectories[0],
                    init_func=init)
    plt.show()


def animation_2(xx, nt, trajectories,legends=None,colors=['k-','b-']):
    
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln1, = ax.plot([], [], colors[0], animated=True)
    ln2, = ax.plot([], [], colors[1], animated=True)
    
    def init():
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)
        if legends!=None:
            ax.legend(legends)
        return ln1,

    def update(i):
        xdata=xx 
        ydata=trajectories[0][i]
        ln1.set_data(xdata,ydata)
        ydata=trajectories[1][i]
        ln2.set_data(xdata,ydata)
        return ln1,ln2,

    ani = FuncAnimation(fig, update, np.arange(nt),
                            init_func=init, blit=True)
    plt.show()


def animation_3(xx, nt, trajectories,legends=None,colors=['k-','b-','r-']):

    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln1, = ax.plot([], [], colors[0], animated=True)
    ln2, = ax.plot([], [], colors[1], animated=True)
    ln3, = ax.plot([], [], colors[2], animated=True)
    
    def init():
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)
        if legends!=None:
            ax.legend(legends)
        return ln1,

    def update(i):
        xdata=xx 
        ydata=trajectories[0][i]
        ln1.set_data(xdata,ydata)
        ydata=trajectories[1][i]
        ln2.set_data(xdata,ydata)
        ydata=trajectories[2][i]
        ln3.set_data(xdata,ydata)
        return ln1,ln2,ln3,

    ani = FuncAnimation(fig, update, np.arange(nt),
                            init_func=init, blit=True)
    plt.show()


def animation_4(xx, nt, trajectories,legends=None,colors=['k-','b-','r-','g-']):

    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln1, = ax.plot([], [], colors[0], animated=True)
    ln2, = ax.plot([], [], colors[1], animated=True)
    ln3, = ax.plot([], [], colors[2], animated=True)
    ln4, = ax.plot([], [], colors[3], animated=True)
    
    def init():
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)
        if legends!=None:
            ax.legend(legends)
        return ln1,

    def update(i):
        xdata=xx 
        ydata=trajectories[0][i]
        ln1.set_data(xdata,ydata)
        ydata=trajectories[1][i]
        ln2.set_data(xdata,ydata)
        ydata=trajectories[2][i]
        ln3.set_data(xdata,ydata)
        ydata=trajectories[3][i]
        ln4.set_data(xdata,ydata)
        return ln1,ln2,ln3,ln4,

    ani = FuncAnimation(fig, update, np.arange(nt),
                            init_func=init, blit=True)
    plt.show()
