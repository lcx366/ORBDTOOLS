import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares

from ...utils.math import getangle 
from ..radar.gibbs import gibbs_assist,gibbs,herrick_gibbs
from ..common import slant,get_a0

def EHnu(nu,ecc,degrees=True):
    """
    Calculate (E)ccentric anomaly or (H)yperbolic anomaly from the true anomaly ν.

    Usage:
        >>> cosE,sinE = EHnu(50,0.1)
        >>> coshH,sinhH = EHnu(50,1.1)
    Inputs:
        nu -> [float] True anomaly, [deg] or [rad]
        ecc -> [float] Eccentricity
        degrees -> [bool,optional,default=True] Unit of true anomaly. If True, ν is in [deg]; otherwose in [rad].
    Outputs:
        (cosE,sinE) for ecc < 1 or (coshH,sinhH) for ecc > 1    
    """
    if degrees: nu = np.deg2rad(nu)
    d = 1 + ecc*np.cos(nu)
    if ecc < 1: # Eccentric anomaly
        cosE = (ecc + np.cos(nu)) /  d
        sinE = np.sqrt(1-ecc**2)*np.sin(nu)/d
        return cosE,sinE  
    else: # Hyperbolic anomaly
        coshH = (ecc + np.cos(nu)) /  d
        sinhH = np.sqrt(ecc**2-1)*np.sin(nu)/d
        return coshH,sinhH          

def doubleR_assist(r12,mu,xyz_site3p,los3p,tof,tau,degrees=True):
    """
    Auxiliary function for least squares estimation of r12 in Double-R method.

    Usage:
        >>> F1,F2,ele = doubleR_assist(r12,mu,xyz_site3p,los3p,tof,tau)
    Inputs:
        r12 -> [list of float with two elements] Position vector of space objects at the start and end time moment of the observation in GCRF
        mu -> [float] GM of the central body of attraction
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        F1 -> [float] O-C of t1-t2
        F2 -> [float] O-C of t3-t2
        ele -> [list] Classical orbital elements with values as follows
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    """
    tau1,tau2 = -tau[0],tau[1]

    R1_vec,R2_vec,R3_vec = xyz_site3p
    los1,los2,los3 = los3p

    rho12,r12_vec = slant(r12,los3p[:2],xyz_site3p[:2])
    r1_vec,r2_vec = r12_vec

    W_vec = np.cross(r1_vec,r2_vec)
    rho3 = np.dot(-R3_vec,W_vec)/np.dot(los3,W_vec)
    r3_vec = R3_vec + los3*rho3
    r123_vec = np.vstack([r12_vec,r3_vec])

    method,mode,nu21,nu32,nu31 = gibbs_assist(r123_vec,degrees)

    if method == 'Gibbs':
        ele = gibbs(mu,r123_vec,degrees)
    elif method == 'Herrick-Gibbs':
        ele = herrick_gibbs(mu,r123_vec,tof,tau,degrees) 

    a,ecc,inc,raan,argp,nu2 = ele    

    nu1 = nu2 - nu21 
    nu3 = nu2 + nu32

    if ecc < 1: # for elliptical orbit

        cosE1,sinE1  = EHnu(nu1,ecc,degrees)
        cosE2,sinE2  = EHnu(nu2,ecc,degrees)
        cosE3,sinE3  = EHnu(nu3,ecc,degrees)

        sinE32 = sinE3*cosE2 - cosE3*sinE2
        cosE32 = cosE3*cosE2 + sinE3*sinE2
        sinE21 = sinE2*cosE1 - cosE2*sinE1
        cosE21 = cosE2*cosE1 + sinE2*sinE1

        E32 = getangle(sinE32,cosE32) # in radian
        E21 = getangle(sinE21,cosE21) # in radian

        M32 = E32 - ecc*(sinE3 - sinE2)
        M12 = -E21 - ecc*(sinE1 - sinE2)

        n = np.sqrt(mu/a**3)

    else: # for hyperbolic orbit

        coshH1,sinhH1  = EHnu(nu1,ecc)
        coshH2,sinhH2  = EHnu(nu2,ecc)
        coshH3,sinhH3  = EHnu(nu3,ecc)

        sinhH32 = sinhH3*coshH2 - coshH3*sinhH2
        coshH32 = coshH3*coshH2 - sinhH3*sinhH2
        sinhH21 = sinhH2*coshH1 - coshH2*sinhH1
        coshH21 = coshH2*coshH1 - sinhH2*sinhH1

        H32 = np.arcsinh(sinhH32)
        H21 = np.arcsinh(sinhH21)

        M32 = -H32 + ecc*(sinhH3 - sinhH2)
        M12 = H21 + ecc*(sinhH1 - sinhH2)

        n = np.sqrt(-mu/a**3)

    F1 = tau1 - M12/n
    F2 = tau2 - M32/n

    return F1,F2,ele

def fun_resi(r12,mu,xyz_site3p,los3p,tof,tau,degrees):
    """
    Residual function for least squares estimation of r12 in Double-R method.

    Usage:
        >>> residuals = fun_resi(r12,xyz_site3p,los3p,tof,tau,degrees)
    Inputs:
        r12 -> [list of float with two elements] Position vector of space objects at the start and end time moment of the observation in GCRF
        mu -> [float] GM of the central body of attraction
        xyz_site3p -> [2D array with shape of 3x3] Cartesian coordinates of the site at three moments in GCRF
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        degrees -> [bool] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        residuals -> Norm of O-C of (t1-t2,t3-t2)
    """
    F1,F2,ele = doubleR_assist(r12,mu,xyz_site3p,los3p,tof,tau,degrees)
    residuals = norm([F1,F2])

    return residuals

def doubleR_estimate(mu,tof,tau,los3p,xyz_site3p,degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Double-R method.

    Usage:
        >>> ele = doubleR_estimate(mu,tof,tau,los3p,xyz_site3p)
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at three moments
        xyz_site3p_nd -> [2D array with shape of nx3] Cartesian coordinates of the site at three moments
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
    """
    # Guess the initial semi-major axis
    a0 = get_a0(mu,los3p,xyz_site3p,tof)
    r12 = np.array([a0]*2)
    res = least_squares(fun_resi,r12, args=(mu,xyz_site3p,los3p,tof,tau,degrees),loss='huber',bounds=([1,1],[8,8]),method='dogbox')
    F1,F2,ele = doubleR_assist(res.x,mu,xyz_site3p,los3p,tof,tau,degrees)
     
    return ele   