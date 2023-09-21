import numpy as np
from numpy.linalg import norm

from ...utils import Const
from ...transform.kep_rv_trans import rv2coe

def gibbs(mu,r123_vec,degrees=True):
    """
    Compute the classical orbital elements at Median epoch from three position vectors using Gibbs Method.

    Usage:
        >>> # Example 7-3. Using Three Position Vectors for GIBBS Orbit Determination in <Fundamentals of Astrodynamics and Applications (Fourth Edition)> by David A. Vallado.
        >>> import numpy as np
        >>> r123_vec = np.array([[0,0,6378.137],[0,-4464.696,-5102.509],[0,5740.323,3189.068]]) # Position vectors of space object in [km]
        >>> mu = 398600.4418 # GM, [km^3/s^2]  
        >>> from orbdtools.iod.radar.gibbs import gibbs
        >>> ele = gibbs(mu,r123_vec)
        >>> print(ele)
        >>> # [6.63981040e+03 4.08085652e-02 9.00000000e+01 9.00000000e+01 1.05664222e+02 1.23149848e+02]
    Inputs:
        mu -> [float] GM of the central body of attraction
        r123_vec -> [2D array of float] Three position vectors in form of [[x1,y1,z1],...,[x3,y3,z3]]
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [1D array with 6 elements] Classical orbital elements with
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    Note: The Gibbs Method is superior for seperation of position vectors above five degrees.  
    """
    r1_vec,r2_vec,r3_vec = r123_vec
    r1,r2,r3  = r123 = norm(r123_vec,axis=1)
    r1_uec,r2_uec,r3_uec = r123_vec/r123[:,None]

    c12 = np.cross(r1_vec, r2_vec)
    c23 = np.cross(r2_vec, r3_vec)
    c31 = np.cross(r3_vec, r1_vec)
    c23u = c23/norm(c23)

    # Checking colplanarity
    alpha = 90 - np.degrees(np.arccos(np.dot(c23u,r1_uec)))
    if alpha > 1: 
        raise Exception('Coplanar condition for vectors is not satisfied!') 
        
    N_vec = r1*c23 + r2*c31 + r3*c12   
    N = norm(N_vec)

    D_vec = c12 + c23 + c31
    D = norm(D_vec)
        
    S_vec = r1_vec*(r2 - r3) + r2_vec*(r3 - r1) + r3_vec*(r1 - r2)
    v2_vec = np.sqrt(mu/(N*D))*(np.cross(D_vec,r2_uec) + S_vec)

    ele = rv2coe(np.hstack([r2_vec,v2_vec]),mu,degrees)

    return ele

def herrick_gibbs(mu,r123_vec,tof,tau,degrees=True):
    """
    Compute the classical orbital elements at Median epoch from three position vectors using Herrick-Gibbs Method.

    Usage:
        >>> # Example 7-4. Using Three Position Vectors for HERRICK-GIBBS Orbit Determination in <Fundamentals of Astrodynamics and Applications (Fourth Edition)> by David A. Vallado.
        >>> import numpy as np
        >>> from astropy.time import Time
        >>> r123_vec = np.array([[3419.85564,6019.82602,2784.60022],[2935.91195,6326.18324,2660.59584],[2434.95202,6597.38674,2521.52311]]) # Position vectors of the space object
        >>> t123 = Time(['2016-01-02T03:04:05.000','2016-01-02T03:05:21.480','2016-01-02T03:06:38.040']) # Time of the three position vectors in UTC
        >>> tof = (t123[-1] - t123[0]).sec
        >>> tau = [(t123[1] - t123[0]).sec,(t123[2] - t123[1]).sec]
        >>> mu = 398600.4418 # GM, [km^3/s^2]  
        >>> from orbdtools.iod.radar.gibbs import herrick_gibbs
        >>> ele = herrick_gibbs(mu,r123_vec,tof,tau)
        >>> print(ele)
        >>> # [8.29125159e+03 9.99644405e-02 2.49999950e+01 3.00000022e+02 1.18000168e+02 4.49981113e+00]
    Inputs:
        mu -> [float] GM of the central body of attraction
        r123_vec -> [2D array of float] Three position vectors in form of [[x1,y1,z1],...,[x3,y3,z3]]
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [1D array with 6 elements] Classical orbital elements with
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    Note: The Herrick-Gibbs Method is superior for seperation of position vectors below one degree.    
    """
    r1_vec,r2_vec,r3_vec = r123_vec
    r1,r2,r3  = r123 = norm(r123_vec,axis=1)
    r1_uec,r2_uec,r3_uec = r123_vec/r123[:,None]

    c23 = np.cross(r2_vec, r3_vec)
    c23u = c23/norm(c23)

    # Checking colplanarity
    alpha = 90 - np.degrees(np.arccos(np.dot(c23u,r1_uec)))
    if alpha > 1: 
        raise Exception('Coplanar condition for vectors is not satisfied!')  
        
    v2_vec_p1 = -tau[1]*(1/(tau[0]*tof) + mu/(12*r1**3))*r1_vec
    v2_vec_p2 = (tau[1] - tau[0])*(1/(tau[0]*tau[1]) + mu/(12*r2**3))*r2_vec
    v2_vec_p3 = tau[0]*(1/(tau[1]*tof) + mu/(12*r3**3))*r3_vec
    v2_vec = v2_vec_p1 + v2_vec_p2 + v2_vec_p3

    ele = rv2coe(np.hstack([r2_vec,v2_vec]),mu,degrees)

    return ele

def gibbs_assist(r123_vec,degrees=True):
    """
    Judge whether to use the method Gibbs or Herrick-Gibbs to determine the initial orbit based on three position vectors.

    Usage:
        >>> method,mode,nu21,nu32,nu31 = gibbs_assist(r123_vec)
    Inputs:
        r123_vec -> [2D array of float] Three position vectors in form of [[x1,y1,z1],...,[x3,y3,z3]]
    Outputs:
        method -> [str] 'Herrick-Gibbs' or 'Gibbs'
        mode -> [str] Available modes include 'r1-r2r3p', 'r1-r2r3r', 'r1r2r3p', 'r1r2r3r', 'r1r2-r3p', and 'r1r2-r3r'.
            'r1-r2r3p' means prograde orbit with r1 and r2 constitute a reflex angles in time order;
            'r1-r2r3r' means retrograde orbit with r1 and r2 constitute a reflex angles in time order;
            'r1r2r3p' means prograde orbit with r1, r2, and r3 constitute a normal angular sequence in time order;
            'r1r2r3r' means retrograde orbit with r1, r2, and r3 constitute a normal angular sequence in time order;
            'r1r2-r3p' means prograde orbit with r2 and r3 constitute a reflex angles in time order;
            'r1r2-r3r' means retrograde orbit with r2 and r3 constitute a reflex angles in time order;
        nu21 -> [float] Angle swept by the vector from t1 to t2 in [deg]
        nu32 -> [float] Angle swept by the vector from t2 to t3 in [deg]
        nu31 -> [float] Angle swept by the vector from t1 to t3 in [deg]
    """
    twopi = Const.twopi

    r1_vec,r2_vec,r3_vec = r123_vec

    # unit vector
    r1_uec,r2_uec,r3_uec = r123_vec/norm(r123_vec,axis=1)[:,None]

    c12 = np.cross(r1_vec, r2_vec)
    c23 = np.cross(r2_vec, r3_vec)
    c31 = np.cross(r3_vec, r1_vec)

    if np.dot(c23,c31) > 0:
        if c23[2] > 0:
            mode = 'r1-r2r3p' # prograde orbit
        else:
            mode = 'r1-r2r3r' # retrograde orbit
        nu21 = twopi - np.arccos(np.dot(r2_uec,r1_uec)) # nu2 - nu1
        nu32 = np.arccos(np.dot(r3_uec,r2_uec))
        nu31 = twopi - np.arccos(np.dot(r3_uec,r1_uec))      
    elif np.dot(c12,c23) > 0:
        if c12[2] > 0:
            mode = 'r1r2r3p'
        else:
            mode = 'r1r2r3r'  
        nu21 = np.arccos(np.dot(r2_uec,r1_uec))
        nu32 = np.arccos(np.dot(r3_uec,r2_uec))  
        nu31 = np.arccos(np.dot(r3_uec,r1_uec))
    elif np.dot(c31,c12) > 0: 
        if c31[2] > 0:
            mode = 'r1r2-r3p'
        else:
            mode = 'r1r2-r3r' 
        nu21 = np.arccos(np.dot(r2_uec,r1_uec)) 
        nu32 = twopi - np.arccos(np.dot(r3_uec,r2_uec)) 
        nu31 = twopi - np.arccos(np.dot(r3_uec,r1_uec))

    nu_max = np.rad2deg(max(nu21,nu32))
    if degrees: nu21,nu32,nu31 = np.rad2deg([nu21,nu32,nu31])    

    if nu_max < 5: 
        method = 'Herrick-Gibbs'
    else:
        method = 'Gibbs'    
  
    return method,mode,nu21,nu32,nu31            
             
def gibbs_estimate(mu,tof,tau,xyz3p,degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from three-points radar measurements using Gibbs/Herrick-Gibbs Method.

    Usage:
        >>> import numpy as np
        >>> from astropy.time import Time
        >>> r123_vec = np.array([[3419.85564,6019.82602,2784.60022],[2935.91195,6326.18324,2660.59584],[2434.95202,6597.38674,2521.52311]])
        >>> t123 = Time(['2016-01-02T03:04:05.000','2016-01-02T03:05:21.480','2016-01-02T03:06:38.040'])
        >>> tof = (t123[-1] - t123[0]).sec
        >>> tau = [(t123[1] - t123[0]).sec,(t123[2] - t123[1]).sec]
        >>> mu = 398600.4418 # GM, [km^3/s^2]  
        >>> from orbdtools.iod.radar.gibbs import gibbs_estimate
        >>> ele = gibbs_estimate(mu,tof,tau,r123_vec)
        >>> print(ele)
        >>> # [8.29125159e+03 9.99644405e-02 2.49999950e+01 3.00000022e+02 1.18000168e+02 4.49981113e+00]
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        xyz3p -> [2D array of float] Position vectors of the space object in form of [[x1,y1,z1],...,[x3,y3,z3]]
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [1D array with 6 elements] Classical orbital elements with
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    """
    method,mode,nu21,nu32,nu31 = gibbs_assist(xyz3p,degrees)
    
    if method == 'Gibbs':
        ele = gibbs(mu,xyz3p,degrees)
    elif method == 'Herrick-Gibbs':
        ele = herrick_gibbs(mu,xyz3p,tof,tau,degrees) 
    return ele