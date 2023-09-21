from numpy.linalg import det
from numpy.polynomial import Polynomial
import numpy as np

from ...transform.kep_rv_trans import rv2coe

def laplace_iod_root(mu,t_,xyz_sitenp,losnp):
    """
    Find the real, positive roots of the 8-th order polynomial equation.

    Usage:
        >>> poly_roots_real_positive,params = laplace_iod_root(mu,t_,xyz_sitenp,losnp)
    Inputs:
        mu -> [float] GM of the central body of attraction
        t_ -> [array of float] Elapsed time since the Intermediate Epoch
        xyz_sitenp -> [2D array with shape of nx3] Cartesian coordinates of the site
        losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
    Outputs:
        poly_roots_real_positive -> [array of float] Real, positive roots for the 8-th order polinominal equation, which is also the distance of the space object w.r.t. the geocenter at the Intermediate epoch of the observation
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
    """
    t2 = 0
    
    # For sapce objects
    # Fit a 2-order polynomial
    order = 2 # empirical value
    losx = Polynomial.fit(t_, losnp[:,0],order).convert()
    losy = Polynomial.fit(t_, losnp[:,1],order).convert()
    losz = Polynomial.fit(t_, losnp[:,2],order).convert()

    losx_d1 = losx.deriv(m=1)
    losy_d1 = losy.deriv(m=1)
    losz_d1 = losz.deriv(m=1)

    losx_d2 = losx.deriv(m=2)
    losy_d2 = losy.deriv(m=2)
    losz_d2 = losz.deriv(m=2)

    los2 = np.array([losx(t2),losy(t2),losz(t2)])
    los2_d1 = np.array([losx_d1(t2),losy_d1(t2),losz_d1(t2)])
    los2_d2 = np.array([losx_d2(t2),losy_d2(t2),losz_d2(t2)])

    # For sites
    # Fit a 4-order polynomial
    obsx = Polynomial.fit(t_, xyz_sitenp[:,0], order).convert()
    obsy = Polynomial.fit(t_, xyz_sitenp[:,1], order).convert()
    obsz = Polynomial.fit(t_, xyz_sitenp[:,2], order).convert()

    obsx_d1 = obsx.deriv(m=1)
    obsy_d1 = obsy.deriv(m=1)
    obsz_d1 = obsz.deriv(m=1)

    obsx_d2 = obsx.deriv(m=2)
    obsy_d2 = obsy.deriv(m=2)
    obsz_d2 = obsz.deriv(m=2)

    R2_vec = np.array([obsx(t2),obsy(t2),obsz(t2)])
    R2_vec_d1 = np.array([obsx_d1(t2),obsy_d1(t2),obsz_d1(t2)])
    R2_vec_d2 = np.array([obsx_d2(t2),obsy_d2(t2),obsz_d2(t2)])
        
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
    
    # get the real, positive roots of the 8-th order polynomial equation
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
    Compute the state vector of the space object at the Intermediate epoch of the observation from optical angle-only measurement data using Multiple-points Laplace method.

    Usage:
        >>> rv2_vec = laplace_iod(mu,params,r2)
    Inputs:
        mu -> [float] GM of the central body of attraction 
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
        r2 -> [float] Distance of the space object w.r.t. the geocenter at the Intermediate epoch of the observation, which is solved from the 8-th order polinominal equation 
    Outputs:
        rv2_vec -> [array of float] State vector(position and velocity) of the space object at the Intermediate epoch of the observation
    """
    D,D1,D2,D3,D4,los2,los2_d1,R2_vec,R2_vec_d1 = params
    rho2 = -2*(D1/D + mu/r2**3*D2/D) # rho2 must be greater than 0
    rho2_d1 = -(D3/D + mu/r2**3*D4/D)
    r2_vec = R2_vec + rho2*los2
    v2_vec = rho2_d1*los2 + rho2*los2_d1 + R2_vec_d1
    rv2_vec = np.hstack([r2_vec,v2_vec])
    
    return rv2_vec

def multilaplace_estimate(mu,t_,xyz_sitenp,losnp,degrees=True):
    """
    Estimate the classical orbital elements at Intermediate epoch from optical angle-only measurements using Multiple-points Laplace method.

    Usage:
        >>> eles = multilaplace_estimate(mu,t_,xyz_sitenp,losnp)
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
        Bate R R, Mueller D D, White J E, et al. Fundamentals of astrodynamics(2nd)[M]. Courier Dover Publications, 2020.     
    """
    poly_roots_real_positive,params = laplace_iod_root(mu,t_,xyz_sitenp,losnp)
     
    eles = []         
    for r2 in poly_roots_real_positive:  
        rv2_vec = laplace_iod(mu,params,r2)
        ele = rv2coe(rv2_vec,mu,degrees)
        eles.append(ele)
    return eles