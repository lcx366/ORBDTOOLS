import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares

from ...transform.kep_rv_trans import coe2rv     
from ..lambert.izzo import izzo_iod  
from ..lambert.universal import universal_iod
from ..common import slant,get_a0,coe_propagation 

def gooding_assist(rho13,mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees):
    """
    Auxiliary function for least squares estimation of rho13 in Gooding method.

    Usage:
        >>> rho2_vec,eles_t2 = gooding_assist(rho13,mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees)
    Inputs:
        mu -> [float] GM of the central body of attraction
        rho13 -> [list of float with two elements] Slant distance of space objects w.r.t. site at the start and end time moment of the observation
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        tm -> [int,optional,default=1] Lambert transfer mode. If tm = 1, then short way transfer mode; else if tm = -1, long way transfer mode.  
        M -> [int,optional,default=0] Number of full revolutions in transfer
        method -> [str,optional,default='universal'] Method for solving the Lambert's problem. Available options include 'universal' and 'izzo'
        degrees -> [bool] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        rho2_vec -> Slant Vector at the Median epoch
        eles_t2 -> [list] Classical orbital elements at the Median epoch with values as follows
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    """

    rho1,rho3 = rho13
    R1_vec,R2_vec,R3_vec = xyz_site3p
    los1,los2,los3 = los3p
    r1_vec = R1_vec + los1*rho1
    r3_vec = R3_vec + los3*rho3

    if method == 'universal':
        eles_t1 = universal_iod(mu,r1_vec,r3_vec,tof,tm,degrees)
    elif method == 'izzo':    
        eles_t1 = izzo_iod(mu,r1_vec,r3_vec,tof,tm,M,degrees)[0]
    elif method == 'battin':
        pass

    # propagate the orbital elements
    eles_t2 = coe_propagation(eles_t1,tau[0],mu,degrees)
    r2_vec = coe2rv(eles_t2,mu,degrees)[:3]   
    rho2_vec = r2_vec - R2_vec

    return rho2_vec,eles_t2

def fun_resi(rho13,mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees):
    """
    Residual function for least squares estimation of rho13 in Gooding method.

    Usage:
        >>> residuals = fun_resi(rho13,mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees)
    Inputs:
        rho13 -> [list of float with two elements] Slant distance of space objects w.r.t. site at the start and end time moment of the observation
        mu -> [float] GM of the central body of attraction
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        tm -> [int,optional,default=1] Lambert transfer mode. If tm = 1, then short way transfer mode; else if tm = -1, long way transfer mode.  
        M -> [int,optional,default=0] Number of full revolutions in transfer
        method -> [str,optional,default='universal'] Method for solving the Lambert's problem. Available options include 'universal' and 'izzo'
        degrees -> [bool] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        residuals -> Difference between one and cosine of the angle spaned by the observational LOS and the calculated LOS at the Median epoch
    """
    rho2_vec,eles_t2 = gooding_assist(rho13,mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees)
    proj_insight = np.dot(rho2_vec/norm(rho2_vec),los3p[1])
    residuals = proj_insight - 1

    return residuals    

def gooding_estimate(mu,tof,tau,los3p,xyz_site3p,tm=1,M=0,method='universal',degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from optical angle-only measurement data using Gooding method.

    Usage:
        >>> ele = gooding_estimate(mu,tof,tau,los3p,xyz_site3p)
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        xyz_site3p -> [2D array with shape of nx3] Cartesian coordinates of the site at three moments
        tm -> [int,optional,default=1] Lambert transfer mode. If tm = 1, then short way transfer mode; else if tm = -1, long way transfer mode.  
        M -> [int,optional,default=0] Number of full revolutions in transfer
        method -> [str,optional,default='universal'] Method for solving the Lambert's problem. Available options include 'universal' and 'izzo'
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [list] Classical orbital elements with values as follows
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    References: 
        Gooding R H. A new procedure for the solution of the classical problem of minimal orbit determination from three lines of sight[J]. Celestial Mechanics and Dynamical Astronomy, 1996, 66: 387-423.      
        Izzo D. Revisiting Lambert’s problem[J]. Celestial Mechanics and Dynamical Astronomy, 2015, 121: 1-15.
        Curtis H D. Orbital Mechanics for Engineering Students: Revised 4th edition[M]. Butterworth-Heinemann, 2020.     
    """
    # Guess the initial semi-major axis
    a0 = get_a0(mu,los3p,xyz_site3p,tof)
    r13 = np.array([a0]*2)
    rho13,r13_vec = slant(r13,los3p[::2],xyz_site3p[::2]) 

    res = least_squares(fun_resi,rho13, args=(mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees),loss='huber',bounds=([0,0],[8,8]),method='dogbox')
    rho2_vec,ele = gooding_assist(res.x,mu,xyz_site3p,los3p,tof,tau,tm,M,method,degrees)

    return ele  