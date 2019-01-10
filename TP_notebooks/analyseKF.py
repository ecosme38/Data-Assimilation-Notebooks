import scipy.linalg as lin
import numpy as np

def analyseKF(up,Sp,H,yo,R):
    # Kalman filter analysis
  
    #Kalman gain
    HS = np.dot(H,Sp)
    Kg = np.dot( np.dot(Sp,HS.T) ,
                 lin.inv( np.dot(HS,HS.T) + R ) )
    # Analysis
    uu =  up + np.dot(Kg,
                      yo-np.dot(H,up) )
    P  =  np.dot( Sp-np.dot(Kg,HS) ,
                  Sp.T )
    P  = 0.5 * (P+P.T)          # we force symmetry
    S  = lin.sqrtm(P)         # square root decomposition of Pf

    return uu, S
