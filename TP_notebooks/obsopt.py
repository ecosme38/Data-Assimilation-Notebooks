import numpy as np

class Obsopt:

    def __init__(self,nx,xsub,nt,tsub):
        self.nx = nx
        self.nt = nt
        self.tsub = tsub # time and
        self.xsub = xsub # space subsampling
        self.yo={}       # observation vectors
        
        obsxmask = np.array([i%xsub==0 for i in range(1,nx+1)],dtype=float)

        obs_location = np.flatnonzero(obsxmask)
        self.nobs   = np.size(obs_location)
        
        self.mat=np.zeros((self.nobs,nx))
        for i in range(self.nobs):
            self.mat[i,obs_location[i]]=1.
        

    def isobserved(self,t):
        return t%self.tsub==0
        
    def gen_obs(self,model,u0,sigmao):
        true=[u0] # true trajectory
        u=u0
        for t in range(self.nt):
            if self.isobserved(t):
                noise = np.random.normal(0.,sigmao,u.size)
                self.yo[t] = np.dot(self.mat,(u + noise))
            u = model.step(u)
            true.append(u)

        if self.isobserved(self.nt):
            noise = np.random.normal(0.,sigmao,u.size)
            self.yo[self.nt] = np.dot(self.mat,(u + noise))

        return true
    
    def dir(self,t,u):
        if self.isobserved(t):
            return np.dot(self.mat,u)

    def tan(self,t,u):
        if self.isobserved(t):
            return np.dot(self.mat,u)

    def adj(self,t,y):
        if self.isobserved(t):
            return np.dot(self.mat.T,y)

    def misfit(self,t,u):
         if self.isobserved(t):
             return np.dot(self.mat,u) - self.yo[t]
       
