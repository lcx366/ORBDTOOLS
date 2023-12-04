import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares
from scipy.linalg import lstsq

from ...transform.kep_rv_trans import rv2coe

def D_M(f,g,xyz_sitenp,losnp):
    """
    Design matrix for equation A * rv_vec = B

    Usage:
        >>> A,B = D_M(f,g,xyz_sitenp,losnp)
    Inputs:
        f -> [array_like] Lagrangian coefficient f for r_vec = f(t,r0_vec,v0_vec)*r0_vec + g(t,r0_vec,v0_vec)*v0_vec
        g -> [array_like] Lagrangian coefficient g for r_vec = f(t,r0_vec,v0_vec)*r0_vec + g(t,r0_vec,v0_vec)*v0_vec
        xyz_sitenp -> [2D array with shape of nx3] Cartesian coordinates of the site at three moments
        losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
    Outputs:
        A -> [array with shape of nx3x6] Design matrix
        B -> [array with shape of nx3x1] B in equation of A * rv_vec = B
    """
    
    lx,ly,lz = losnp.T
    xs,ys,zs = xyz_sitenp.T

    zeros = np.zeros_like(f)

    A = np.array([[f*lz,zeros,-f*lx,g*lz,zeros,-g*lx],\
        [zeros,f*lz,-f*ly,zeros,g*lz,-g*ly],\
        [f*ly,-f*lx,zeros,g*ly,-g*lx,zeros]])
    B = np.array([[lz*xs-lx*zs,lz*ys-ly*zs,ly*xs-lx*ys]])  

    A,B = A.transpose(2,0,1),B.T 

    return A,B  

def fun_resi(rv_vec,mu,xyz_sitenp,losnp,tau):
    """
    Residual function for least squares estimation of rv_vec in FG-Series method.

    Usage:
        >>> residuals = fun_resi(rv_vec,mu,xyz_sitenp,losnp,tau)
    Inputs:
        rv_vec -> [list of float with two elements] State vector(position and velocity) of the space object at the intermediate moment of the observation
        mu -> [float] GM of the central body of attraction
        xyz_sitenp -> [2D array with shape of nx3] Cartesian coordinates of the site
        losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
        tau -> [array of float] Elapsed time since the Intermediate Epoch, namely t-t_ref, where t_ref takes the intermediate moment of the observation
    Outputs:
        residuals -> Norm of O-C of A * rv_vec - B
    """
    r_vec,v_vec = rv_vec[:3],rv_vec[3:]
    r = norm(r_vec)

    u = mu/r**3
    p = np.dot(r_vec,v_vec)/r**2
    q = np.dot(v_vec,v_vec)/r**2 - u

    f = 1 - u*tau**2/2 + u*p*tau**3/2 + u*(u - 15*p**2 + 3*q)*tau**4/24 + u*p*(7*p**2 - u - 3*q)*tau**5/8 
    g = tau - u*tau**3/6 + u*p*tau**4/4 + u*(u - 45*p**2 + 9*q)*tau**5/120

    A,B = D_M(f,g,xyz_sitenp,losnp)
    residuals = norm(A@rv_vec - B[:,:,-1],axis=1)

    return residuals      

def fg_series_optical(mu,t_,xyz_sitenp,losnp,degrees=True,rv0=None):
    """
    Estimate the classical orbital elements at Intermediate epoch from optical angle-only measurements using FG-Series method.

    Usage:
        >>> ele = fg_series_optical(mu,t_,xyz_sitenp,losnp)
    Inputs:
        mu -> [float] GM of the central body of attraction
        t_ -> [array of float] Elapsed time since the Intermediate Epoch
        xyz_sitenp -> [2D array with shape of nx3] Cartesian coordinates of the site
        losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [list] Classical orbital elements with values as follows
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    references: 
        李光宇.天体测量和天体力学基础[M].科学出版社,2015.   
        刘林.卫星轨道力学算法[M].南京大学出版社,2019.
        Bate R R, Mueller D D, White J E, et al. Fundamentals of astrodynamics(2nd)[M]. Courier Dover Publications, 2020.     
    """

    if rv0 is None:
        f0,g0 = np.ones_like(t_),t_ # Initial Lagrangian Coefficient f and g
        A,B = D_M(f0,g0,xyz_sitenp,losnp)
        num_eqs = len(t_)*3
        rv0,_resi,_rnk,_s = lstsq(A.reshape(num_eqs,6),B.reshape(num_eqs))
    res = least_squares(fun_resi,rv0, args=(mu,xyz_sitenp,losnp,t_),method='lm')
    ele = rv2coe(res.x,mu,degrees)

    return ele