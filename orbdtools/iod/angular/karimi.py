import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares
from scipy.linalg import lstsq

from ...utils import Const
from ..radar.multifitting import ellipse_fitting

# Design matrix
def D_M(c1,c3,pos_obs,los):
    '''
    Considering the number of parameters to be fitted, this method is only suitable for a small number of data points. 
    For a large number of data points, please use the method of liulin.
    f and g -> [array_like] Lagrangian coefficient
    los -> [nx3 array] : Line-Of-Sight(LOS) vectors for targre w.r.t. site in ICRF at n moments

    reference: Karimi R R, Mortari D. Initial orbit determination using multiple observations[J]. Celestial Mechanics and Dynamical Astronomy, 2011, 109(2): 167-180.
    '''

    n = len(los)
    A_e,B_e = [],[]
    for i in range(n-2):
        A = np.zeros((3,n))
        los1,los2,los3 = los[i:i+3]
        A[:,i],A[:,i+1],A[:,i+2] = c1[i]*los1,-los2,c3[i]*los3
        A_e.append(A)

        B = pos_obs[i+1] - c1[i]*pos_obs[i] - c3[i]*pos_obs[i+2]
        B_e.append(B)
        
    A_e = np.vstack(A_e) 
    B_e = np.hstack(B_e) 

    return A_e,B_e 

def fun_resi(rho,pos_obs_nd,los,tau1,tau3):

    mu = Const.mu_nd # GM for center of attraction

    r_vec = pos_obs_nd + rho[:,None]*los
    u = mu/norm(r_vec[1:-1])**3

    f1 = 1 - u*tau1**2/2
    f3 = 1 - u*tau3**2/2
    g1 = tau1 - u*tau1**3/6
    g3 = tau3 - u*tau3**3/6
    fg = f1*g3 - f3*g1
    c1,c3 = g3/fg,-g1/fg

    A,B = D_M(c1,c3,pos_obs_nd,los)

    print(rho)

    residuals = A@rho - B

    return residuals        

def karimi_estimate(tnp,pos_obs_nd,losnp,degrees=True):

    t_ = (tnp - tnp[0]).sec/Const.T_nd
    t_diff = np.diff(t_)
    tau1,tau3 = np.array([-t_diff[:-1],t_diff[1:]])

    f1 = f3 = 1
    g1,g3 = tau1,tau3
    fg = f1*g3 - f3*g1
    c1,c3 = g3/fg,-g1/fg

    A,B = D_M(c1,c3,pos_obs_nd,losnp)

    rho0,_resi,_rnk,_s = lstsq(A,B)
    res = least_squares(fun_resi,rho0, args=(pos_obs_nd,losnp,tau1,tau3))
    r_vec = pos_obs_nd + res.x[:,None]*losnp

    tn,ele = ellipse_fitting(tnp,r_vec,degrees)

    return tn,ele