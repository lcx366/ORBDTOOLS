import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares
from scipy.linalg import lstsq
from astropy.coordinates import cartesian_to_spherical

from ...transform.kep_rv_trans import rv2coe
from ...utils import Const
from ..common import slant,get_a0

# Design matrix
def D_M(f,g,alpha_vec,delta_vec,r_vec):
    '''
    f and g -> [array_like] Lagrangian coefficient
    los -> [nx3 array] : Line-Of-Sight(LOS) vectors for targre w.r.t. site in ICRF at n moments
    '''

    A1 = np.hstack([f[:,None]*alpha_vec,g[:,None]*alpha_vec])
    A2 = np.hstack([f[:,None]*delta_vec,g[:,None]*delta_vec])
    B1 = np.sum(r_vec * alpha_vec,axis=1)  
    B2 = np.sum(r_vec * delta_vec,axis=1) 

    A = np.vstack([A1,A2]) 
    B = np.hstack([B1,B2])

    return A,B

def fun_resi(rv,R_vec,losnp,tau):

    mu = Const.mu_nd # GM for center of attraction

    r_vec,v_vec = rv[:3],rv[3:]
    r = norm(r_vec)

    u = mu/r**3
    p = np.dot(r_vec,v_vec)/r**2
    q = np.dot(v_vec,v_vec)/r**2 - u

    f = 1 - u*tau**2/2 + u*p*tau**3/2 + u*(u - 15*p**2 + 3*q)*tau**4/24 + u*p*(7*p**2 - u - 3*q)*tau**5/8 
    g = tau - u*tau**3/6 + u*p*tau**4/4 + u*(u - 45*p**2 + 9*q)*tau**5/120

    r_vecs = f[:,None]*r_vec + g[:,None]*v_vec
    rho_vecs = r_vecs - R_vec
    rhos = norm(rho_vecs,axis = 1)

    los_ref = rho_vecs/rhos[:,None]

    _,dec_ref,ra_ref = cartesian_to_spherical(los_ref[:,0],los_ref[:,1],los_ref[:,2])
    alpha_vec = np.array([-np.sin(ra_ref),np.cos(ra_ref),np.zeros_like(f)]).T
    delta_vec = np.cross(los_ref,alpha_vec)

    A,B = D_M(f,g,alpha_vec,delta_vec,r_vecs)

    residuals = A@rv - B

    return residuals      

def ref_vec_estimate(pos_obs3p,los3p,R_vec,losnp,tof,tnp,alpha_vec,delta_vec,degrees=True):

    '''
    reference : 贾沛璋.初轨计算的参考矢量法[J].天文学报,1997.
    '''

    mu = Const.mu_nd

    tau = ((tnp - tnp[-1]) + (tnp - tnp[0])).sec/2/Const.T_nd

    f0,g0 = np.ones_like(tau),tau
    # compute the major axis of orbit
    a = get_a0(los3p,pos_obs3p,tof)
    rho,r_vec = slant(a,losnp,R_vec)

    A,B = D_M(f0,g0,alpha_vec,delta_vec,r_vec)

    rv0,_resi,_rnk,_s = lstsq(A,B)

    res = least_squares(fun_resi,rv0,args=(R_vec,losnp,tau),loss='huber')
    ele = rv2coe(res.x,mu,degrees)

    return ele