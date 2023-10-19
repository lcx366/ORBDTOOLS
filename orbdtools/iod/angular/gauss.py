import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares

from ...transform.kep_rv_trans import rv2coe

def gauss_iod_root(mu,tof,tau,xyz_site3p,los3p):
    """
    Find the real, positive roots of the 8-th order polynomial equation.

    Usage:
        >>> poly_roots_real_positive,params = gauss_iod_root(mu,tof,tau,xyz_site3p,los3p)
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
    Outputs:
        poly_roots_real_positive -> [array of float] Real, positive roots for the 8-th order polinominal equation, which is also the distance of the space object w.r.t. the geocenter at the Median epoch
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
    """
    los1,los2,los3 = los3p
    R1_vec,R2_vec,R3_vec = xyz_site3p

    tau1,tau3 =  -tau[0],tau[1]
    tau = tof

    p = np.array((np.zeros((3,)),np.zeros((3,)),np.zeros((3,))))
    p1 = np.cross(los2, los3)
    p2 = np.cross(los1, los3)
    p3 = np.cross(los1, los2)
    p = np.stack([p1,p2,p3]).T

    D0  = np.dot(los1, p1)
    D = np.dot(xyz_site3p,p)

    A = (-D[0,1]*tau3/tau + D[1,1] + D[2,1]*tau1/tof)/D0
    B = (D[0,1]*(tau3**2-tof**2)*tau3/tof + D[2,1]*(tof**2-tau1**2)*tau1/tof)/(6*D0)

    E = np.dot(R2_vec, los2)
    R2_square = np.dot(R2_vec, R2_vec)

    a = -(A**2 + 2*A*E + R2_square)
    b = -2*mu*B*(A+E)
    c = -(mu*B)**2

    # get the real, positive roots of the 8-th order polynomial equation
    poly_coeffs = np.zeros(9)
    poly_coeffs[0] = 1
    poly_coeffs[2] = a
    poly_coeffs[5] = b
    poly_coeffs[8] = c

    poly_roots = np.roots(poly_coeffs)
    poly_roots_real_positive = np.real(poly_roots[np.isreal(poly_roots)&(poly_roots > 0)])

    params = A,B,D,D0,tau1,tau3,tof
    
    return poly_roots_real_positive,params   
        
def gauss_iod(mu,xyz_site3p,los3p,params,r2):
    """
    Compute the preliminary state vector of space object at Median epoch from optical angle-only measurements using Gauss method.

    Usage:
        >>> rv2_vec = gauss_iod(mu,xyz_site3p,los3p,params,r2)
    Inputs:
        mu -> [float] GM of the central body of attraction
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
        r2 -> [float] Distance of the space object w.r.t. the geocenter at the Median epoch, which is solved from the 8-th order polinominal equation 
    Outputs:
        rv2_vec -> [array of float] State vector(position and velocity) of the space object at the Median epoch for 3-points observation
    """
    R1_vec,R2_vec,R3_vec = xyz_site3p
    los1,los2,los3 = los3p
    A,B,D,D0,tau1,tau3,tof = params

    rho2 = A + mu*B/r2**3
    
    rho1_p1 = 6*(D[2,0]*tau1/tau3 + D[1,0]*tof/tau3)*r2**3 + mu*D[2,0]*(tof**2-tau1**2)*tau1/tau3
    rho1_p2 = 6*r2**3 + mu*(tof**2 - tau3**2)
    rho1 = (rho1_p1/rho1_p2 - D[0,0])/D0

    rho3_p1 = 6*(D[0,2]*tau3/tau1 - D[1,2]*tof/tau1)*r2**3 + mu*D[0,2]*(tof**2-tau3**2)*tau3/tau1
    rho3_p2 = 6*r2**3 + mu*(tof**2 - tau1**2)
    rho3 = (rho3_p1/rho3_p2 - D[2,2])/D0

    rhos = rho1,rho2,rho3 

    r1_vec = R1_vec + los1*rho1
    r2_vec = R2_vec + los2*rho2
    r3_vec = R3_vec + los3*rho3

    u = mu/r2**3
    f1 = 1 - u*tau1**2/2
    f3 = 1 - u*tau3**2/2
    g1 = tau1 - u*tau1**3/6
    g3 = tau3 - u*tau3**3/6

    v2_vec = (-f3*r1_vec + f1*r3_vec)/(f1*g3 - f3*g1)
    rv2_vec = np.hstack([r2_vec,v2_vec])

    return rv2_vec

def fun_resi(rv2_vec,mu,xyz_site3p,los3p,params):
    """
    Residual function for least squares estimation of rv2_vec in Gauss method.

    Usage:
        >>> residuals = fun_resi(rv2_vec,mu,xyz_site3p,los3p,params)
    Inputs:
        rv2_vec -> [array of float] State vector(position and velocity) of the space object at the Median epoch for 3-points observation
        mu -> [float] GM of the central body
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        params -> [list of float] Parameters associated to the 8-th order polinominal equation
    Outputs:
        residuals -> Norm of O-C of rv2_vec
    """
    R1_vec,R2_vec,R3_vec = xyz_site3p
    los1,los2,los3 = los3p
    A,B,D,D0,tau1,tau3,tof = params

    r2_vec,v2_vec = rv2_vec[:3],rv2_vec[3:]
    r2 = norm(r2_vec)
    u = mu/r2**3
    p = np.dot(r2_vec,v2_vec)/r2**2
    q = np.dot(v2_vec,v2_vec)/r2**2 - u

    f1 = 1 - u*tau1**2/2 + u*p*tau1**3/2 + u*(u - 15*p**2 + 3*q)*tau1**4/24 + u*p*(7*p**2 - u - 3*q)*tau1**5/8
    f3 = 1 - u*tau3**2/2 + u*p*tau3**3/2 + u*(u - 15*p**2 + 3*q)*tau3**4/24 + u*p*(7*p**2 - u - 3*q)*tau3**5/8
    g1 = tau1 - u*tau1**3/6 + u*p*tau1**4/4 + u*(u - 45*p**2 + 9*q)*tau1**5/120
    g3 = tau3 - u*tau3**3/6 + u*p*tau3**4/4 + u*(u - 45*p**2 + 9*q)*tau3**5/120

    fg = f1*g3 - f3*g1
    c1,c3 = g3/fg,-g1/fg

    rho1 = (-D[0,0] + D[1,0]/c1 - c3*D[2,0]/c1)/D0
    rho2 = (-c1*D[0,1] + D[1,1] - c3*D[2,1])/D0
    rho3 = (-c1/c3*D[0,2] + D[1,2]/c3 - D[2,2])/D0

    r1_vec = R1_vec + los1*rho1
    r2_vec = R2_vec + los2*rho2
    r3_vec = R3_vec + los3*rho3

    v2_vec = (-f3*r1_vec + f1*r3_vec)/fg
    
    residuals = rv2_vec - np.hstack([r2_vec,v2_vec])

    return residuals   

def gauss_estimate(mu,tof,tau,los3p,xyz_site3p,degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Gauss method.

    Usage:
        >>> eles = gauss_estimate(mu,tof,tau,los3p,xyz_site3p)
    Inputs:
        mu -> [float] GM of the central body of attraction
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
        Curtis H D. Orbital Mechanics for Engineering Students: Revised 4th edition[M]. Butterworth-Heinemann, 2020.           
    """
    poly_roots_real_positive,params = gauss_iod_root(mu,tof,tau,xyz_site3p,los3p)

    eles = []
    for r2 in poly_roots_real_positive:
        rv2_vec0 = gauss_iod(mu,xyz_site3p,los3p,params,r2)
        res = least_squares(fun_resi,rv2_vec0, args=(mu,xyz_site3p,los3p,params),method='lm')
        ele = rv2coe(res.x,mu,degrees) 
        eles.append(ele)       

    return eles