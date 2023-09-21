import numpy as np
from numpy.linalg import norm
from scipy.optimize import fsolve

from ...transform.kep_rv_trans import rv2coe

def stumpff(z):
    """
    The Stumpff functions c2(z),c3(z) developed by Karl Stumpff.
    """
    if z > 0:
        c2 = (1 - np.cos(np.sqrt(z))) / z
        c3 = (np.sqrt(z) - np.sin(np.sqrt(z))) / z**1.5
    elif z < 0:
        c2 = -(np.cosh(np.sqrt(-z) - 1)) / z
        c3 = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (-z)**1.5
    else:
        c2 = 0.5
        c3 = 1/6

    return c2,c3

def func_lambert(z,*args):
    """
    Find z from the Lambert equation.
    """
    mu,r1_plus_r2,A,tof = args

    c2,c3 = stumpff(z)
    y = r1_plus_r2 + A * (z * c3 - 1) / np.sqrt(c2)
    xi = np.sqrt(y / c2)
    return (xi ** 3 * c3 + A * np.sqrt(y)) - np.sqrt(mu)*tof

def universal_iod(mu,r1_vec,r2_vec,tof,tm=1,degrees=True):
    """
    Solve the Lambert's problem using the universal variable approach.

    Inputs:
        mu -> [float] GM of the center of attraction
        r1_vec -> The beginning position vector of space object in unit of [L_nd]
        r2_vec -> The ending position vector of object in unit of [L_nd]
        tof -> Time of flight in unit of [T_nd]
        tm -> [int,optional,default=1] Transfer mode. If tm = 1,then short way transfer mode; else if tm = -1, long way transfer mode.  
        degrees -> [bool,optional,default=True] Units for angle variables in orbital elements estimated. If True, angle variables in [deg], otherwise in [rad] 
    Outputs:
        ele -> [list] Classical orbital elements with values as follows
            a -> [float] Semi-major axis in [L_nd]
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    Note:
        Multiple revolutions are not supported.
    """

    r1,r2 = norm([r1_vec,r2_vec],axis=1)
    r1_times_r2 = r1 * r2
    r1_plus_r2 = r1 + r2

    cos_dnu = np.dot(r1_vec,r2_vec) / r1_times_r2

    A = tm * np.sqrt(r1_times_r2 * (1 + cos_dnu))
    if A == 0: raise Exception("Cannot compute orbit, phase angle is 180 degrees")

    args = (mu,r1_plus_r2,A,tof)
    z0 = 1 # initial guess of z
    x = fsolve(func_lambert, z0, args)
    z = x[0]   

    c2,c3 = stumpff(z)
    y = r1_plus_r2 + A * (z * c3 - 1) / np.sqrt(c2)

    # calculate the Lagrange functions
    f = 1 - y / r1
    g = A * np.sqrt(y / mu)
    gdot = 1 - y / r2
 
    # Calculate pair of velocity solutions
    v1_vec = (r2_vec - f * r1_vec) / g
    v2_vec = (gdot * r2_vec - r1_vec) / g

    ele = rv2coe(np.hstack([r1_vec,v1_vec]),mu,degrees) 

    return ele