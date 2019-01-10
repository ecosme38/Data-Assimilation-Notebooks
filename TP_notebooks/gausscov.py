import scipy.linalg as lin
from math import exp
import numpy as np

class gausscov:

    def __init__(self,nx,sigma,L, indic):

        self.mat = np.zeros((nx,nx))

        for j in range(nx):
            for i in range(nx):
                # matB(i,j)=sigma*sigma*exp(-(((i-j)*dx)**2)/(2.*L*L));   	
                k=min(abs(j-i),nx-j+i)
                dx=1./nx
                self.mat[i,j]=sigma*sigma*exp(-((k*dx)**2)/(2.*L*L))

        self.mat=0.5*(self.mat+self.mat.T) # ensure symetry

        self.factor(indic)
        
    def factor(self,indic) :
        if indic==1 :
            self.inv=lin.inv(self.mat)
        elif indic==2 :
            self.sqr=lin.sqrtm(self.mat)
        elif indic==3 :
            self.ext=cholesky(self.mat)
        else: 
            raise ValueError('unknown indic in gausscov')
    
