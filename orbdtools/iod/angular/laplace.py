from scipy.interpolate import KroghInterpolator
from numpy.linalg import det,norm
import numpy as np

from ...transform.kep_rv_trans import rv2coe

def laplace_iod_root(mu,tof,tau,xyz_site3p,los3p):
    """
    Find the real, positive roots of the 8-th order polynomial equation.

    Usage:
        >>> poly_roots_real_positive,params = laplace_iod_root(mu,tof,tau,xyz_site3p,los3p)
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
    Outputs:
        poly_roots_real_positive -> [array of float] Real, positive roots for the 8-th order polinominal equation, which is also the distance of the space object w.r.t. the geocenter at the Median moment
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
    """
    t_ = t1,t2,t3 = [0,tau[0],tof]
    los1,los2,los3 = los3p
    R1_vec,R2_vec,R3_vec = xyz_site3p
    
    # for space object
    poly_tar = KroghInterpolator(t_,los3p)
    los2_d1,los2_d2 = poly_tar.derivatives(t2,3)[1:]

    # for site
    poly_obs = KroghInterpolator(t_,xyz_site3p)
    R2_vec_d1,R2_vec_d2 = poly_obs.derivatives(t2,3)[1:]
        
    D1 = det([los2,los2_d1,R2_vec_d2])
    D2 = det([los2,los2_d1,R2_vec])
    D3 = det([los2,R2_vec_d2,los2_d2])
    D4 = det([los2,R2_vec,los2_d2])

    D = 2*det([los2,los2_d1,los2_d2])
    A = -2*D1/D
    B = -2*mu*D2/D
    C = np.dot(los2,R2_vec)  
    
    a = -(A**2+2*A*C+np.dot(R2_vec,R2_vec))
    b = -2*B*(A+C)
    c = -B**2
    
    # get real, positive roots of the 8-th order polynomial equation
    poly_coeffs = np.zeros(9)
    poly_coeffs[0] = 1
    poly_coeffs[2] = a
    poly_coeffs[5] = b
    poly_coeffs[8] = c

    poly_roots = np.roots(poly_coeffs)
    poly_roots_real_positive = np.real(poly_roots[np.isreal(poly_roots)&(poly_roots > 0)])
    
    params = D,D1,D2,D3,D4,los2,los2_d1,R2_vec,R2_vec_d1
    
    return poly_roots_real_positive,params
        
def laplace_iod(mu,params,r2):
    """
    Compute the state vector of the space object at Median epoch from optical angle-only measurement data using Laplace method.

    Usage:
        >>> rv2_vec = laplace_iod(mu,params,r2)
    Inputs:
        mu -> [float] GM of the central body of attraction
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
        r2 -> [float] Distance of the space object w.r.t. the geocenter at the Median epoch, which is solved from the 8-th order polinominal equation 
    Outputs:
        rv2_vec -> [array of float] State vector(position and velocity) of the space object at the Median epoch for 3-points observation
    """
    D,D1,D2,D3,D4,los2,los2_d1,R2_vec,R2_vec_d1 = params
    rho2 = -2*(D1/D + mu/r2**3*D2/D) # rho2 must be greater than 0
    rho2_d1 = -(D3/D + mu/r2**3*D4/D)
    r2_vec = R2_vec + rho2*los2
    v2_vec = rho2_d1*los2 + rho2*los2_d1 + R2_vec_d1
    rv2_vec = np.hstack([r2_vec,v2_vec])
    
    return rv2_vec  

def laplace_estimate(mu,tof,tau,los3p,xyz_site3p,degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Laplace method.

    Usage:
        >>> eles = laplace_estimate(mu,tof,tau,los3p,xyz_site3p)
    Inputs:
        mu -> [float] GM of the center of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        xyz_site3p -> [2D array with shape of nx3] Cartesian coordinates of the site at three moments
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
        Vallado D. Fundamentals of Astrodynamics and Applications(4th)[M], Microcosm Press, 2013. 
        Bate R R, Mueller D D, White J E, et al. Fundamentals of astrodynamics(2nd)[M]. Courier Dover Publications, 2020.            
    """
    poly_roots_real_positive,params = laplace_iod_root(mu,tof,tau,xyz_site3p,los3p)
             
    eles = []
    for r2 in poly_roots_real_positive:  
        rv2_vec = laplace_iod(mu,params,r2)
        ele = rv2coe(rv2_vec,mu,degrees)
        eles.append(ele)     
    return eles         