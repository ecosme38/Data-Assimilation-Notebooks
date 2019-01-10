import numpy as np
from scipy.linalg import inv

class Variational:
    def __init__(self,ubkg=None, nt=None, B=None, M=None, H=None, R=None, precond=True):
        self.prec=precond
        self.B=B
        self.M=M
        self.H=H
        self.Rinv=inv(R)
        self.ubkg=ubkg
        self.nt = nt

    def cost(self,v):

        f = self.simvar(v,1)

        return f
        
    def grad(self,v):

        g = self.simvar(v,2)

        return g
        
    def simvar(self,v,indic):
        
        # Change of variable if precond
        if self.prec :
            u  = self.B.sqr.dot(v) + self.ubkg
            gb = v        # gradient of background term
            Jb = v.dot(v) # cost of background term
        else:
            u  = v + self.ubkg
            gb = np.dot(self.B.inv,v) # gradient of background term
            Jb = np.dot(v,gb)         # cost of background term

        u_trj=list()                    # Storage of reference trajectory

        # Time Loop. Cost function evaluation
        Jo=0.
        for it in range(self.nt):
            u_trj.append(u)
            if self.H.isobserved(it):
                misfit=self.H.misfit(it,u) # d=Hx-xobs
                Jo=Jo+misfit.dot(self.Rinv.dot(misfit))
            u=self.M.step(u)

        if self.H.isobserved(self.nt):
            misfit=self.H.misfit(self.nt,u) # d=Hx-xobsx
            Jo=Jo+misfit.dot(self.Rinv.dot(misfit))
            
        J=0.5*(Jb+Jo) # Total cost function
        # print 'J: ',J
        if indic==2 :
            # reverse time loop, Gradient evaluation

            uad=np.zeros(self.M.nx)

            if self.H.isobserved(self.nt):
                misfit=self.H.misfit(self.nt,u) # d=Hx-xobs
                uad = uad + self.H.adj(self.nt,self.Rinv.dot(misfit))

            for itr in reversed(range(self.nt)):
                # Retreive chekpoints
                u  = u_trj[itr]
                # One backward step
                uad=self.M.step_adj(u,uad);
                # Calculation of adjoint forcing
                if self.H.isobserved(itr):
                    misfit=self.H.misfit(itr,u) # d=Hx-xobs
                    uad = uad + self.H.adj(itr,self.Rinv.dot(misfit))

            # Adjoint of the change of varable, if needed 
            if self.prec :
                g=np.transpose(self.B.sqr).dot(uad) + gb # total gradient
            else:
                g=uad + gb
            # print 'G: ',g.dot(g)
            return g

        else:

            return J
