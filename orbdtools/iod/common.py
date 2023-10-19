import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar,least_squares
import pandas as pd

from ..utils import Const
from ..transform.kep_rv_trans import coe2rv
from ..transform.orbele_trans import nu_to_Me,nu_to_Mh,Me_to_nu,Mh_to_nu
        
def coe_propagation(coe0,t_,mu,degrees):
    """
    Propagation of classical orbital elements with the two-body unperturbed dynamic model.

    Usage:
        >>> coe_t = coe_propagation(coe0,t_,mu,degrees)
    Inputs:
        coe0 -> [list of float with 6 elements] Classical orbital elements at time of t0
        t_ -> [array of float] Elapsed time since t0
        mu -> [float] GM of the central body of attraction
        degrees -> [bool] Unit of angular variables in classical orbital elements. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        coe_t -> [list of float with 6 elements] Classical orbital elements at the time of elapsed
    """
    twopi = Const.twopi

    a0,ecc0,inc0,raan0,argp0,nu0 = coe0 
    n = np.sqrt(mu/np.abs(a0)**3) # mean motion

    a_t = np.ones_like(t_)*a0
    ecc_t = np.ones_like(t_)*ecc0
    inc_t = np.ones_like(t_)*inc0
    raan_t = np.ones_like(t_)*raan0
    argp_t = np.ones_like(t_)*argp0

    if ecc0 < 1: # For elliptical orbit
        M0 = nu_to_Me(nu0,ecc0,degrees)

        if degrees:
            M_t = M0 + np.rad2deg(n*t_)
            nu_t = Me_to_nu(M_t,ecc0,True)%360
        else:
            M_t = M0 + n*t_    
            nu_t = Me_to_nu(M_t,ecc0,False)%twopi

    else: # For hyperbolic orbit
        M0 = nu_to_Mh(nu0,ecc0,degrees)    

        if degrees:
            M_t = (M0 + np.rad2deg(n*t_))%360
            nu_t = Mh_to_nu(M_t,ecc0,True)%360
        else:
            M_t = (M0 + n*t_)%twopi    
            nu_t = Mh_to_nu(M_t,ecc0,False)%twopi

    coe_t = np.stack([a_t,ecc_t,inc_t,raan_t,argp_t,nu_t]).T     

    return coe_t

def func_circular(a,*args):
    """
    Auxiliary function for least squares estimation of semi-major axis in Near-Circular Orbit Hypothesis.

    Usage:
        >>> resi = func_circular(a,*args)
    Inputs:
        a -> [float] Semi-major axis
        args -> optional arguments include
            mu -> [float] GM of the central body of attraction
            los_vec -> [array_like] Unit vector of the Line-Of-Sight(LOS) from the site to the space object
            R_vec -> [array_like] Position vector of the site
            tof -> [float] Time of flight
    Outputs:
        resi -> Magnitude of O-C for dot(r1_vec,r2_vec)
    """
    mu,los_vec,R_vec,tof = args
    rho,r_vec = slant(a,los_vec,R_vec)
    n = np.sqrt(mu/a**3)
    resi = np.abs(np.dot(r_vec[0],r_vec[-1]) - a**2*np.cos(n*tof)) 
    
    return resi    

def slant(r,los_vec,R_vec):
    """
    Calculate the slant distance of the space object w.r.t. the site.

    usage:
        >>> rho,r_vec = slant(r,los_vec,R_vec)
    Inputs:
        r -> [array_like] Distance to the central attractor
        los_vec -> [array_like] Unit vector of the Line-Of-Sight(LOS) from the site to the space object
        R_vec -> [array_like] Position vector of the site
    Outputs:
        rho -> [array_like] Slant distance of the space object w.r.t. the site
        r_vec -> [array_like] Position vector of the space object
    """
    if R_vec.ndim == 1:
        R = norm(R_vec)
        C = 2*np.dot(los_vec,R_vec)
        q = C**2 - 4*(R**2 - r**2)
        if q < 0: q = 0
        rho = (-C + np.sqrt(q))/2     
        r_vec = rho*los_vec + R_vec
    else:    
        R = norm(R_vec,axis=1)
        C = 2*(los_vec * R_vec).sum(axis=1)
        q = C**2 - 4*(R**2 - r**2)
        q[q < 0] = 0
        rho = (-C + np.sqrt(q))/2 
        r_vec = rho[:,None]*los_vec + R_vec
    return rho,r_vec

def get_a0(mu,los_vec,R_vec,tof):
    """
    Guess the initial radius of the circular orbit.

    Usage:
        >>> a0 = get_a0(mu,los_vec,R_vec,tof)
    Inputs:
        mu -> [float] GM of the central body of attraction
        los_vec -> [array_like] Unit vector of the Line-Of-Sight(LOS) from the site to the space object
        R_vec -> [array_like] Position vector of the site
        tof -> [float] Time of flight
    Outputs:
        a0 -> [float] Initial radius of the circular orbit    
    """
    C = np.sum(los_vec*R_vec,axis=1)
    if C.max() < 0: 
        r_min = norm(R_vec - los_vec*C[:,None],axis=1).max()
    else:
        r_min = norm(R_vec,axis=1).max()  

    args = (mu,los_vec,R_vec,tof)
    res = minimize_scalar(func_circular, bounds=(r_min,r_min*1.5), args = args,method='bounded')

    if not res.success: raise Exception(res.message)           

    return res.x    

def fun_resi_optical(coe0,t_,mu,losnp,xyz_sitenp,degrees):
    """
    Calculate the residual of IOD for optical angle-only measurement data.

    Usage:
        >>> residuals = fun_resi_optical(coe0,t_,mu,losnp,xyz_sitenp,degrees)
    Inputs:
        coe0 -> [list of float with 6 elements] Classical orbital elements at the reference epoch
        t_ -> [array of float] Elapsed time since the reference epoch
        mu -> [float] GM of the central body of attraction
        losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
        xyz_sitenp -> [2D array with shape of nx3] Cartesian coordinates of the site
        degrees -> [bool] Unit of angular variables in classical orbital elements. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        residuals -> Difference between one and cosine of the angle spaned by the observed LOS and the calculated LOS
    """
    coe_t = coe_propagation(coe0,t_,mu,degrees)
    rs_vec = coe2rv(coe_t,mu,degrees)[:,:3]
    rhos_vec = rs_vec - xyz_sitenp
    proj_insight = np.sum(rhos_vec*losnp,axis=1)
    residuals = proj_insight/norm(rhos_vec,axis=1) - 1

    return residuals            

def fun_resi_radar(coe0,t_,mu,xyz_np,degrees):
    """
    Calculate the residual of IOD for radar range+angle measurements.

    Usage:
        >>> residuals = fun_resi_radar(coe0,t_,mu,xyz_np,degrees)
    Inputs:
        coe0 -> [list of float with 6 elements] Classical orbital elements at the reference epoch
        mu -> [float] GM of the central body of attraction
        t_ -> [array of float] Elapsed time since the reference epoch
        xyz_np -> [2D array with shape of nx3] Cartesian coordinates of the space object
        degrees -> [bool] Unit of angular variables in classical orbital elements. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        residuals -> [Array of float] Norm of the O-C of the position vector of the space object
    """

    coe_t = coe_propagation(coe0,t_,mu,degrees)
    rs_vec = coe2rv(coe_t,mu,degrees)[:,:3]

    residuals = norm(rs_vec - xyz_np,axis=1)

    return residuals

def getnpoints(n,t,*args):
    """
    Extract n approximately uniformly distributed data points. 

    Usage:
        >>> t_np,xyz_site_np,radec_np = getnpoints(t,xyz_site,radec)
    Inputs:
        t -> [Astropy Time or Skyfield Time] Time sequence of observation
        xyz_site -> [array like, optional] Cartesian coordinates of the site
        radec -> [array like, optional] RA and Dec of the space object
    Outputs:
        t_np -> [Astropy Time or Skyfield Time] Time extracted
        xyz_site_np -> [array like, optional] Cartesian coordinates extracted
        radec_np -> [array like, optional] RA and Dec extracted
    """
    if n > 0:
        tof = t[-1] - t[0]
        indices = np.zeros(n,dtype=int)

        for i in range(n):
            mt = t[0] + i*tof/(n-1) # time moment of the node of arc
            mi = np.argmin(np.abs(t - mt)) # index corresponding to the node of arc
            indices[i] = mi  

        indices = np.unique(indices) # remove duplicate indices
        out = [t[indices]]       

        for cont in args:
            out.append(cont[indices])    
    else:
        # for all data
        out = [t]   
        for cont in args:
            out.append(cont)  

    return out 

def to_ele_dict_optical(mu,t0,t_,eles,losnp,xyz_sitenp,degrees,rms_tol=1e-8):   
    """
    Determine the validity of the orbital elements according to the O-C for optical angle-only measurement data.

    Usage:
        >>> ele_dict = to_ele_dict_optical(mu,t0,t_,ele,losnp,xyz_sitenp,degrees)
    Inputs:
        mu -> [float] GM of central body of attraction
        t0 -> [Astropy Time or Skyfield Time] Reference epoch of the orbital elements
        t_ -> [array of float] Elapsed time since the reference epoch t0
        eles -> [list of float with 6 elements] Classical orbital elements at time of t0
        losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
        xyz_sitenp -> [2D array with shape of nx3] Cartesian coordinates of the site
        degrees -> [bool] Unit of angular variables in classical orbital elements. If True, angular variables are in [deg], otherwise in [rad].
        rms_tol -> [float,optional,default=1e-8] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid. 
    Outputs:
        ele_df -> [Dataframe of Pandas] Dataframe for classical orbital elements with keys and values as follows
            epoch -> [str] Epoch of orbital elements in UTC
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
            M -> [float] Mean anomaly, [rad] or [deg]
            h -> [float] Modulus of angular momentum
            status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
            The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        rms -> [float] RMS of O-C     
    """ 
    ele_df = pd.DataFrame({})
    for ele in eles:
        a,ecc,inc,raan,argp,nu = ele 

        if ecc > 1:
            M = nu_to_Mh(nu,ecc,degrees)
        else:    
            M = nu_to_Me(nu,ecc,degrees)

        h = np.sqrt(a * (1 - ecc**2)) 

        # check if the result is valid
        status = 'failed'

        O_C = fun_resi_optical(ele,t_,mu,losnp,xyz_sitenp,degrees)
        rms = np.sqrt(np.dot(O_C,O_C)/len(O_C))
        if rms < rms_tol: status = 'success'
        ele_dict = {'epoch':[t0.isot],'a':[a],'ecc':[ecc],'inc':[inc],'raan':[raan],'argp':[argp],'nu':[nu],'M':[M],'h':[h],'status':[status]}
        ele_df_i = pd.DataFrame(ele_dict)
        ele_df = pd.concat([ele_df, ele_df_i],ignore_index=True)

    return ele_df,rms   

def to_ele_dict_radar(mu,t0,t_,ele,xyz_np,degrees,rms_tol=1e-8):   
    """
    Determine the validity of the orbital elements according to the O-C for radar range+angle measurements.

    Usage:
        >>> ele_dict = to_ele_dict_radar(mu,t0,t_,ele,xyz_np,degrees)
    Inputs:
        mu -> [float] GM of central body of attraction
        t0 -> [Astropy Time or Skyfield Time] Reference epoch of the orbital elements
        t_ -> [array of float] Elapsed time since the reference epoch t0
        eles -> [list of float with 6 elements] Classical orbital elements at time of t0
        xyz_np -> [2D array with shape of nx3] Cartesian coordinates of the space object
        degrees -> [bool] Unit of angular variables in classical orbital elements. If True, angular variables are in [deg], otherwise in [rad].
        rms_tol -> [float,optional,default=1e-8] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid. 
    Outputs:
        ele_df -> [Dataframe of Pandas] Dataframe for classical orbital elements with keys and values as follows
            epoch -> [str] Epoch of orbital elements in UTC
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
            M -> [float] Mean anomaly, [rad] or [deg]
            h -> [float] Modulus of angular momentum
            status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
            The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        rms -> [float] RMS of O-C    
    """ 
    a,ecc,inc,raan,argp,nu = ele 

    if ecc > 1:
        M = nu_to_Mh(nu,ecc,degrees)
    else:    
        M = nu_to_Me(nu,ecc,degrees)

    h = np.sqrt(a * (1 - ecc**2)) 

    # check if the result is valid
    status = 'failed'

    O_C = fun_resi_radar(ele,t_,mu,xyz_np,degrees)
    rms = np.sqrt(np.dot(O_C,O_C)/len(O_C))
    if rms < rms_tol: status = 'success'
    ele_dict = {'epoch':[t0.isot],'a':[a],'ecc':[ecc],'inc':[inc],'raan':[raan],'argp':[argp],'nu':[nu],'M':[M],'h':[h],'status':[status]}
    ele_df = pd.DataFrame(ele_dict)

    return ele_df,rms     