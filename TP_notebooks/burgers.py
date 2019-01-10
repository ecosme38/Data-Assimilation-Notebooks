import numpy as np

class Burgers:
    
    def __init__(self,nx,dx,dt,ns):
        ''' 
        Burgers 1D model
        Entries:
        nx : number of grid points
        dx : space step
        dt : time step
        ns : integer defining the integration scheme:
             0 : Lax-Friedrich 
        '''
        self.nx=nx
        self.ns=ns
        self.cfl=dt/dx
    
    def step(self,up):
        ''' 
        Burgers 1D model
         Entries:
         up : input field
        '''

        # shift up upwind [um1] and downwind [up1] for integration
        up1=np.roll(up,-1)
        um1=np.roll(up,1)
        # u^2/2,x term is centre-discretized:
        B=0.25*self.cfl*(um1*um1-up1*up1)
        if self.ns==0 : # Lax-Friedrich
            un=0.5*(um1+up1)+B
        else:
            raise ValueError('integration scheme?')

        del up1, um1
        return un

    def step_tan(self,up,uptl):
        ''' 
        Burgers 1D tangent linear model
         Entries:
         up : reference direct field
         uptl : input tangent field
        '''

        # shift up upwind [um1] and downwind [up1] for integration
        up1=np.roll(up,-1)
        um1=np.roll(up,1)
        up1tl=np.roll(uptl,-1)
        um1tl=np.roll(uptl,1)
        # u^2/2,x term is centre-discretized:
        Btl=0.5*self.cfl*(umtl1*um1-uptl1*up1)
        if self.ns==0 : # Lax-Friedrich
            uptl=0.5*(umtl1+uptl1)+Btl
        else:
            raise ValueError('integration scheme?')

        del up1tl, um1tl, up1, um1

        return uptl        

    def step_adj(self,up,upad):
        ''' 
        Burgers 1D adjoint model; adjoint of the discrete model
         Entries:
         up : reference direct field
         upad : input adjoint field
        '''
        # shift up upwind [um1] and downwind [up1] for integration
        up1=np.roll(up,-1)
        um1=np.roll(up,1)
        # u^2/2,x term is centre-discretized:
        if self.ns==0 : # Lax-Friedrich
            um1ad=0.5*upad
            up1ad=0.5*upad
            Bad=upad
            upad=np.zeros(self.nx)
        else:
            raise ValueError('integration scheme?')
        um1ad=um1ad+0.5*self.cfl*(um1*Bad)
        up1ad=up1ad-0.5*self.cfl*(up1*Bad)
        Bad=0.
        # shift up upwind (um1) and downwind (up1) for integration
        upad = np.roll(up1ad,1) + np.roll(um1ad,-1)
        
        del up1ad, um1ad, up1, um1

        return upad

    def step_adj_cont(self,up,upad):
        ''' 
        Burgers 1D adjoint model; adjoint of the continuous equation
         Entries:
         up : reference direct field
         upad : input adjoint field
        '''
        # shift up upwind (um1) and downwind (up1) for integration
        up1ad=np.roll(upad,-1)
        um1ad=np.roll(upad,1)
        # u^2/2,x term is centre-discretized:
        B=0.25*self.cfl*up*(up1ad-um1ad);
        if self.ns==0 : # Lax-Friedrich
            unad=0.5*(um1ad+up1ad)-B
        else:
            raise ValueError('integration scheme?')

        del up1ad, um1ad

        return unad
