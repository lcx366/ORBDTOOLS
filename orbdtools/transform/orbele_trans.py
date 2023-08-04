import numpy as np
from scipy.optimize import root_scalar

from ..utils import Const

def _kepler_equation_e(E,M,ecc):
    """
    Kepler Equation and its 1st derivative with respect to Eccentric Anomaly.

    Inputs:
        E -> [array-like,float] Eccentric Anomaly, [radians]
        M -> [array-like,float] True Anomaly, [radians]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
    Outputs:
        Kepler Equation and its derivatives
    """
    res = E_to_M(E,ecc,False) - M, 1 - ecc*np.cos(E)
    return res

def _kepler_equation_h(F,M,ecc):
    """
    Kepler Equation and its 1st derivative with respect to Hyperbolic Anomaly.

    Inputs:
        F -> [array-like,float] Hyperbolic Anomaly, [radians]
        M -> [array-like,float] True Anomaly, [radians]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        Kepler Equation and its derivatives
    """
    res = F_to_M(F,ecc,False) - M, ecc*np.cosh(F) - 1       
    return res

def nu_to_E(nu,ecc,degrees=True):
    """
    Transform to Eccentric Anomaly from True Anomaly.
    
    Usage:
        >>> E = nu_to_E(nu,ecc)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
    Outputs:
        E -> [array-like,float] Eccentric Anomaly, in (-π,π) or (-360,360)
    """
    if degrees:
        nu_rad = np.deg2rad(nu)
        E_rad = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu_rad / 2))
        E = np.rad2deg(E_rad)
    else:
        E = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu / 2))
    return E

def nu_to_F(nu,ecc,degrees=True):
    """
    Transform to Hyperbolic Anomaly from True Anomaly.
    
    Usage:
        >>> F = nu_to_F(nu, ecc)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        F -> [array-like,float] Hyperbolic Anomaly, between -inf and inf, [radians] or [deg]
    """
    if degrees:
        nu_rad = np.deg2rad(nu)
        F_rad = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(nu_rad / 2))
        F = np.rad2deg(E_rad)
    else:    
        F = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(nu / 2))
    return F   

def E_to_M(E,ecc,degrees=True):
    """
    Transform to Mean Anomaly from Eccentric Anomaly.

    Usage:
        >>> M = E_to_M(E, ecc)
    Inputs:
        E -> [array-like,float] Eccentric Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
    Outputs:
        M -> [array-like,float] Mean anomaly, [radians] or [deg]
    """
    if degrees:
        E_rad = np.deg2rad(E)
        M_rad = E_rad - ecc * np.sin(E_rad)
        M = np.rad2deg(M_rad)
    else:
        M = E - ecc * np.sin(E)
    return M  

def F_to_M(F,ecc,degrees=True):
    """
    Transform to Mean Anomaly from Hyperbolic Anomaly.

    Usage:
        >>> M = F_to_M(F, ecc)
    Inputs:
        F -> [array-like,float] Hyperbolic Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        M -> [array-like,float] Mean anomaly, [radians] or [deg]
    """
    if degrees:
        F_rad = np.deg2rad(F)
        M_rad = ecc * np.sinh(F_rad) - F_rad
        M = np.rad2deg(M_rad)
    else:    
        M = ecc * np.sinh(F) - F
    return M

def nu_to_Me(nu,ecc,degrees=True):
    """
    Transform to Mean anomaly from True Anomaly for ellipse trajectories.

    Usage:
        >>> Me = nu_to_Me(nu, ecc)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity
    Outputs:
        Me -> [array-like,float] Mean Anomaly, [radians] or [deg]
    """
    E = nu_to_E(nu,ecc,degrees)
    Me = E_to_M(E,ecc,degrees)

    return Me  

def nu_to_Mp(nu,degrees=True):
    """
    Transform to Mean anomaly from True Anomaly for parabolic trajectories.

    Usage:
        >>> Mp = nu_to_Mp(nu)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc = 1
    Outputs:
        Mp -> [array-like,float] Mean Anomaly, [radians] or [deg]
    """
    if degrees:
        nu_rad = np.deg2rad(nu)
        Mp_rad = np.tan(nu_rad/2)/2 + np.tan(nu_rad/2)**3/6
        Mp = np.rad2deg(Mp_rad)
    else:    
        Mp = np.tan(nu/2)/2 + np.tan(nu/2)**3/6

    return Mp        

def nu_to_Mh(nu,ecc,degrees=True):
    """
    Transform to Mean anomaly from True Anomaly for hyperbolic trahectories.

    Usage:
        >>> Mh = nu_to_Mh(nu, ecc)
    
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        Mh -> [array-like,float] Mean Anomaly, [radians] or [deg]
    """
    F = nu_to_F(nu,ecc,degrees)
    Mh = F_to_M(F,ecc,degrees)

    return Mh      

def M_to_E(M,ecc,degrees=True):
    """
    Transform to Eccentric Anomaly from Mean Anomaly.

    Usage:
        >>> E = M_to_E(M, ecc)
    Inputs:
        M -> [array-like,float] Mean anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1 
    Outputs:
        E -> [array-like,float] Eccentric anomaly, [radians] or [deg]
    Notes:
        This uses a Newton iteration on the Kepler Equation.
    """
    if degrees: 
        M = np.deg2rad(M)
    else:
        M = np.array(M)    
    if M.ndim == 1:
        sols = []
        for M_i,ecc_i in zip(M,ecc):
            sol = root_scalar(_kepler_equation_e,x0=M_i,args=(M_i,ecc_i), fprime=True, method='newton')
            sols.append(sol.root)
        E = np.array(sols)
    elif M.ndim == 0:
        sol = root_scalar(_kepler_equation_e,x0=M,args=(M,ecc), fprime=True, method='newton')
        E = sol.root  
    if degrees: E = np.rad2deg(E)   
    return E  

def M_to_F(M,ecc,degrees=True):
    """
    Transform to Hyperbolic Anomaly from Mean Anomaly.

    Usage:
        >>> F = M_to_F(M, ecc)
    Inputs:
        M -> [array-like,float] Mean anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        F -> [array-like,float] Hyperbolic anomaly
    Notes:
        This uses a Newton iteration on the Kepler Equation.
    """
    if degrees: M = np.deg2rad(M)
    if M.ndim == 1:
        sols = []
        for M_i,ecc_i in zip(M,ecc):
            sol = root_scalar(_kepler_equation_h,x0=M_i,args=(M_i,ecc_i), fprime=True, method='newton')
            sols.append(sol.root)
        F = np.array(sols)
    elif M.ndim == 0:
        sol = root_scalar(_kepler_equation_h,x0=M,args=(M,ecc), fprime=True, method='newton')
        F = sol.root 
    if degrees: F = np.rad2deg(F)   
    return F     

def E_to_nu(E,ecc,degrees=True):
    r"""
    Transform to True anomaly from Eccentric Anomaly.

    Usage:
        >>> nu = E_to_nu(E, ecc)
    Inputs:
        E -> [array-like,float] Eccentric Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1 
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) or (-360,360)
    """
    if degrees:
        E_rad = np.deg2rad(E)
        nu_rad = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E_rad / 2))
        nu = np.rad2deg(nu_rad)
    else:    
        nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2))
    return nu 

def F_to_nu(F, ecc,degrees=True):
    """
    Transform to True anomaly from Hyperbolic Anomaly.

    Usage:
        >>> nu = F_to_nu(F, ecc)
    Inputs:
        F -> [array-like,float] Hyperbolic Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) or (-360,360)
    """
    if degrees:
        F_rad = np.deg2rad(F)
        nu_rad = 2 * np.arctan(np.sqrt((1 + ecc) / (ecc - 1)) * np.tanh(F_rad / 2))
        nu = np.rad2deg(nu_rad)
    else:
        nu = 2 * np.arctan(np.sqrt((1 + ecc) / (ecc - 1)) * np.tanh(F / 2))
    return nu

def Me_to_nu(Me,ecc,degrees=True):
    """
    Transform to True Anomaly from Mean Anomaly for ellipse trajectories.
    
    Usage:
        >>> nu = Me_to_nu(Me, ecc)
    Inputs:
        Me -> [array-like,float] Mean Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
    """
    E = M_to_E(Me,ecc,degrees)
    nu = E_to_nu(E,ecc,degrees) 
    return nu   

def Mp_to_nu(Mp,degrees=True):
    """
    Transform to True Anomaly from Mean Anomaly for parabolic trajectories.
    
    Usage:
        >>> nu = Mp_to_nu(Mp)
    Inputs:
        Mp -> [array-like,float] Mean Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc = 1
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
    """
    if degrees: Mp = np.deg2rad(Mp)
    z = (3*Mp + np.sqrt(1+(3*Mp)**2))**(1/3)
    nu = np.arctan(z-1/z)*2
    if degrees: nu = np.rad2deg(nu)
    return nu        

def Mh_to_nu(Mh,ecc,degrees=True):
    """
    Transform to True Anomaly from Mean Anomaly for hyperbolic trajectories.
    
    Usage:
        >>> nu = Mh_to_nu(Mh, ecc)
    Inputs:
        Mh -> [array-like,float] Mean Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
    """
    F = M_to_F(Mh,ecc,degrees)
    nu = F_to_nu(F,ecc,degrees) 
    return nu        

def coe2nse(a,ecc,inc,raan,argp,nu,degrees=True):
    """
    Convert to non-singular orbital elements from classical orbital elements.
    The non-singular orbital elements exhibit no singularity for zero eccentricity. 

    Usage: 
        >>> a, inc, raan, xi, eta, l = coe2nse(a, ecc, inc, raan, argp, M)
    Inputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
    Outputs:
        a -> [array-like,float] Semi-major axis
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        xi -> [array-like,float] non-singular parameter
        eta -> [array-like,float] non-singular parameter
        l -> [array-like,float] non-singular parameter, [radians] or [deg]  
    """
    M = nu_to_M(nu,ecc,degrees)
    if degrees: 
        argp_rad = np.deg2rad(argp)
        xi = ecc * np.cos(argp_rad)
        eta = ecc * np.sin(argp_rad)
    else:
        xi = ecc * np.cos(argp)
        eta = ecc * np.sin(argp)
    l = argp + M

    return a, inc, raan, xi, eta, l

def nse2coe(a,inc,raan,xi,eta,l,degrees=True):
    """
    Convert to classical orbital elements from non-singular orbital elements.
    The non-singular orbital elements exhibit no singularity for zero eccentricity. 

    Usage:
        >>> a, ecc, inc, raan, argp, M = nse2coe(a, inc, raan, xi, eta, l)
    Inputs:
        a -> [array-like,float] Semi-major axis
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        xi -> [array-like,float] non-singular parameter
        eta -> [array-like,float] non-singular parameter
        l -> [array-like,float] non-singular parameter, [radians] or [deg] 
    Outputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
    """
    ecc = np.sqrt(xi ** 2 + eta ** 2)
    argp = np.arctan2(eta, xi)
    if degrees:
        argp = np.rad2deg(argp) % 360
        M = (l - argp)%360
    else:
        M = (l - argp) % Const.twopi

    nu = M_to_nu(M,ecc,degrees)    
        
    return a, ecc, inc,raan,argp, nu

def coe2mee(a,ecc,inc,raan,argp,nu,degrees=True):
    """
    Convert to modified equinoctial orbital elements from classical orbital elements.
    These direct modified equinoctial equations exhibit no singularity for zero
    eccentricity and orbital inclinations equal to 0 and 90 degrees. However, two of the
    components are singular for an orbital inclination of 180 degrees.

    Usage:
        >>> p, f, g, h, k, L = coe2mee(a, ecc, inc, raan, argp, M)
    Inputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        omega -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
    Outputs:
        p -> [array-like,float] Semi-latus rectum or parameter
        f -> [array-like,float] Equinoctial parameter f
        g -> [array-like,float] Equinoctial parameter g
        h -> [array-like,float] Equinoctial parameter h
        k -> [array-like,float] Equinoctial parameter k
        L -> [array-like,float] Longitude, [radians] or [deg]
    """
    lonper = raan + argp
    p = a*(1-ecc**2)
    M = nu_to_M(nu,ecc,degrees)
    L = lonper + M

    if degrees:
        lonper = np.deg2rad(lonper)
        raan = np.deg2rad(raan)
        inc = np.deg2rad(inc)
    
    f = ecc * np.cos(lonper)
    g = ecc * np.sin(lonper)
    h = np.tan(inc / 2) * np.cos(raan)
    k = np.tan(inc / 2) * np.sin(raan)

    return p, f, g, h, k, L

def mee2coe(p,f,g,h,k,L,degrees=True):
    """
    Convert to classical orbital elements from modified equinoctial orbital elements.

    Usage:
        >>> a, ecc, inc, raan, argp, M = mee2coe(p, f, g, h, k, L)
    Inputs:
        p -> [array-like,float] Semi-latus rectum or parameter
        f -> [array-like,float] Equinoctial parameter f
        g -> [array-like,float] Equinoctial parameter g
        h -> [array-like,float] Equinoctial parameter h
        k -> [array-like,float] Equinoctial parameter k
        L -> [array-like,float] Longitude, [radians] or [deg]
    Outputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        omega -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
    """
    a = p/(1-f**2-g**2)
    ecc = np.sqrt(f ** 2 + g ** 2)
    inc = 2 * np.arctan(np.sqrt(h ** 2 + k ** 2))
    lonper = np.arctan2(g, f)
    raan = np.arctan2(k, h) % Const.twopi
    argp = (lonper - raan) % Const.twopi

    if degrees: 
        lonper = np.rad2deg(lonper)
        inc = np.rad2deg(inc)
        raan = np.rad2deg(raan)
        argp = np.rad2deg(argp)
        M = (L - lonper) % 360
    else:
        M = (L - lonper) % Const.twopi
            
    nu = M_to_nu(M,ecc,degrees)

    return a, ecc, inc, raan, argp, nu 