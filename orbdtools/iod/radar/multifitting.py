import numpy as np
from numpy.linalg import norm,svd
from scipy.integrate import simps
from scipy.optimize import lsq_linear

from ...utils import Const

def ellipse_fitting(mu,tof,t_,r_vec,degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from multiple-points radar measurements using Elliptical Orbit Fitting Method.

    Usage:
        >>> ele = ellipse_fitting(mu,tof,t_,r_vec)
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        t_ -> [array of float] Elapsed time since the Median Epoch
        r_vec -> [2D array of float] Position vectors in form of [[x1,y1,z1],...,[xn,yn,zn]]
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [1D array with 6 elements] Classical orbital elements with
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [radians] or [deg]
            raan -> [float] Longitude of ascending node, [radians] or [deg]
            argp -> [float] Argument of perigee, [radians] or [deg]
            nu -> [float] True anomaly, [radians] or [deg]    
    """
    twopi = Const.twopi

    mi = np.argmin(np.abs(t_)) # index close to the middle moment

    # make an initial guess for h_vec
    h_vec0 = np.mean(np.cross(r_vec[:-1],r_vec[1:]),axis=0)
    h_uvec0 = h_vec0 / norm(h_vec0) # unit vector

    # With the method of Singular Value Decomposition, vt[2] is the normal vector of the best-fit orbit plane
    u, s, vt = svd(r_vec)
    h_uvec = vt[2]/norm(vt[2])
    if np.dot(h_uvec,h_uvec0) < 0: h_uvec = -h_uvec

    # nodel_vec is an unit vector along the line of intersection of the orbit plane and the x-y plane
    nodel_vec = np.cross([0,0,1],h_uvec)

    # p_x and p_y are 2 orthogonal unit vectors on the orbit plane.
    p_x = nodel_vec / norm(nodel_vec)
    p_y = np.cross(h_uvec,p_x)

    # inclination is the angle between the normal vector and the z axis
    inc = np.arccos(np.clip(h_uvec[2],-1,1))

    # raan is the angle between the x axis and the px axis
    raan = np.arctan2(p_x[1],p_x[0])%twopi

    # project all the points onto the plane.
    xs = np.dot(r_vec,p_x)
    ys = np.dot(r_vec,p_y)
    thetas = np.arctan2(ys,xs) # convert to polar coordinates.

    rs = norm(r_vec,axis=1)

    # remove discontinuity
    interupt_flag = np.diff(thetas) < -twopi
    if np.any(interupt_flag): 
        #np.save('thetas.npy',thetas)
        ind, = interupt_flag.nonzero()
        if len(ind) != 1:
            raise Exception('Error in theta series!')
        else:
            thetas[:ind[0]+1] -= twopi

    h_mag = simps(rs**2,thetas)/tof

    p = h_mag**2/mu

    A = np.array([np.cos(thetas),np.sin(thetas)]).T
    res = lsq_linear(A, h_mag**2/mu/rs - 1)
    f,g = res.x

    ecc = np.sqrt(f**2 + g**2)
    a = p/(1 - ecc**2)
    argp = np.arctan2(g,f)%twopi
    nu = (thetas - argp)%twopi

    if degrees:
        inc,raan,argp,nu = np.rad2deg([inc,raan,argp,nu[mi]])

    ele = a,ecc,inc,raan,argp,nu

    return ele
