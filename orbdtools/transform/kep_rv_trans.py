import numpy as np
from numpy.linalg  import norm

from ..utils import Const

def coe2rv(coe,mu,degrees=True):
    """
    Transform classical orbital elements to state vectors.

    Usage:
        >>> import numpy as np
        >>> # classical orbital elements in form of [a, e, i, Ω, ω, ν]
        >>> coe = np.array([7000,0.01,50,100,30,210]) # semi major axis is in [km], and angles are in [deg]
        >>> mu = 398600.4418 # GM for the Reference Earth Model - WGS84, [km^3/s^2] 
        >>> rv = coe2rv(coe,mu)
    Inputs:
        coe -> [array-like] classical orbital elements in form of [a, e, i, Ω, ω, ν], where
            a: semi major axis
            e: eccentricity
            i: inclination, [rad] or [deg]
            Ω: right ascension of the ascending node, [rad] or [deg]
            ω: argument of perigee, [rad] or [deg]
            ν: true anomaly, [rad] or [deg]
        mu -> [float] GM of the central attraction
        degrees -> [bool,optional,default=True] unit of i, Ω, ω, ν
    Outputs:
        rv -> [array-like] state vector in form of [x, y, z, vx, vy, vz] 
    """
    coe_T = np.array(coe).T # transpose 
    if degrees: coe_T[2:] = np.deg2rad(coe_T[2:])  
    rv_T = np.zeros_like(coe_T) 
    a,ecc,inc,raan,argp,nu = coe_T

    p = a * (1 - ecc**2)
    mag_r = p / (1 + ecc * np.cos(nu))

    arglat = argp + nu  # argument of latitude
    sarglat,carglat = np.sin(arglat),np.cos(arglat)

    c4 = np.sqrt(mu / p)
    c5 = ecc * np.cos(argp) + carglat
    c6 = ecc * np.sin(argp) + sarglat

    sinc,cinc = np.sin(inc),np.cos(inc)
    sraan,craan = np.sin(raan),np.cos(raan)

    # position vector
    rv_T[0] = mag_r * (craan * carglat - sraan * cinc * sarglat)
    rv_T[1] = mag_r * (sraan * carglat + cinc * sarglat * craan)
    rv_T[2] = mag_r * sinc * sarglat

    # velocity vector
    rv_T[3] = -c4 * (craan * c6 + sraan * cinc * c5)
    rv_T[4] = -c4 * (sraan * c6 - craan * cinc * c5)
    rv_T[5] = c4 * c5 * sinc
    
    rv = rv_T.T

    return rv

def rv2coe_1d(rv,mu,tol=1e-9):
    """
    Transform state vectors to classical orbital elements.

    Inputs:
        rv -> [1D array of float] state vector in form of [x, y, z, vx, vy, vz]
        mu -> [float] GM of the central attraction
        tol -> [float,optional,default=1e-9] Threshold for small eccentricity or small inclination, where
               if the eccentricity is less than the threshold, it is treated as a circular orbit;
               if the orbital inclination is less than the threshold, it is treated as an equatorial orbit.
    Outputs:
        coe -> [1D array of float] classical orbital elements in form of [a, e, i, Ω, ω, v], where
            a: semi major axis
            e: eccentricity
            i: inclination, [rad]
            Ω: right ascension of the ascending node, [rad]
            ω: argument of perigee, [rad]
            ν: true anomaly, [rad]
    """
    twopi = Const.twopi    

    r,v = rv[:3],rv[3:]   
    coe = np.zeros_like(rv)  

    mag_r,mag_v = norm(r),norm(v)

    h = np.cross(r, v)
    mag_h = norm(h)
    if np.isclose(mag_h,0): return np.empty(6)*np.nan

    n = np.array([-h[1], h[0], 0])
    mag_n = norm(n)
    e = np.cross(v, h)/ mu - r / mag_r
    ecc = mag_e = norm(e)
    p = mag_h**2 / mu
    inc = np.arccos(np.clip(h[2] / mag_h, -1, 1))

    a = 1 / (2 / mag_r - mag_v**2 / mu)

    circular = ecc < tol
    equatorial = np.abs(inc) < tol

    if equatorial and not circular:
        raan = 0
        argp = np.arctan2(e[1], e[0]) % twopi  # Longitude of periapsis
        nu = np.arctan2((h @ np.cross(e, r)) / mag_h, r @ e) % twopi
    elif not equatorial and circular:
        raan = np.arctan2(n[1], n[0]) % twopi
        argp = 0
        # Argument of latitude
        nu = np.arctan2((r @ np.cross(h, n)) / mag_h, r @ n) % twopi
    elif equatorial and circular:
        raan = 0
        argp = 0
        nu = np.arctan2(r[1], r[0]) % twopi  # True longitude
    else:
        nu = np.arccos(np.clip(np.dot(e,r)/(mag_r * mag_e), -1, 1))
        if np.dot(r, v) < 0: 
            nu = twopi - nu
        else:
            nu = nu % twopi    

        raan = np.arccos(np.clip(n[0] / mag_n, -1, 1))
        if n[1] < 0: raan = (twopi - raan) % twopi

        argp = np.arccos(np.clip(np.dot(n, e) / (mag_n * mag_e), -1, 1))
        if e[2] < 0: argp = (twopi - argp) % twopi
        
    coe = np.array([a,ecc,inc,raan,argp,nu])
    
    return coe      

def rv2coe(rvs,mu,degrees=True,tol=1e-9):
    """
    Transform state vectors to classical orbital elements.

    Usage:
        >>> import numpy as np
        >>> rvs = np.array([[ 4.48e+03, -2.79e+03, -4.68e+03,  1.22e+00,6.81e+00, -2.84e+00],[ 5.48e+03, -3.79e+03, -5.68e+03,  1.52e+00,7.81e+00, -3.84e+00]])
        >>> mu = 398600.4418 # GM for the Reference Earth Model - WGS84, [km^3/s^2] 
        >>> coe = rv2coe(rvs,mu)
    Inputs:
        rvs -> [array-like,float] state vector
        mu -> [float] GM of the central attraction
        degrees -> [bool,optional,default=True] unit of i, Ω, ω, v  
        tol -> [float,optional,default=1e-9] Threshold for small eccentricity or small inclination. 
               If the eccentricity is less than the threshold, it is treated as a circular orbit;
               if the orbital inclination is less than the threshold, it is treated as a equatorial orbit.      
    Outputs:
        coe -> [array-like,float] classical orbital elements, where
            a: semi major axis
            e: eccentricity
            i: inclination, [deg] or [rad]
            Ω: right ascension of the ascending node, [deg] or [rad]
            ω: argument of perigee, [deg] or [rad]
            ν: true anomaly, [deg] or [rad]
    """
    rvs = np.array(rvs)

    if rvs.ndim == 1:
        coe = rv2coe_1d(rvs,mu,tol)
    elif rvs.ndim == 2:
        coe = []
        for rv in rvs:
            coe.append(rv2coe_1d(rv,mu,tol))
        coe = np.array(coe) 

    if degrees: 
        coe_T = coe.T
        coe_T[2:] = np.rad2deg(coe_T[2:])
        coe = coe_T.T       
    return coe