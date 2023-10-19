import numpy as np
from numpy.linalg import norm
from scipy.optimize import root_scalar
from scipy.optimize import minimize
from astropy.time import Time
from skyfield.api import EarthSatellite

from .frame_trans import euler2vectors,gcrf_teme_mat
from .kep_rv_trans import rv2coe,coe2rv

from ..utils import Const,data_prepare
from ..utils.math import Matrix_dot_Vector
from ..cod.sgp4.sgp4_init import sgp4init


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
    res = E_to_Me(E,ecc,False) - M, 1 - ecc*np.cos(E)
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
    res = F_to_Mh(F,ecc,False) - M, ecc*np.cosh(F) - 1       
    return res

def nu_to_E(nu,ecc,degrees=True):
    """
    Transform to Eccentric Anomaly from True Anomaly.
    
    Usage:
        >>> E = nu_to_E(nu,ecc)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
        degrees -> [bool,optional,default=True] unit for angular variables 
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
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        F -> [array-like,float] Hyperbolic Anomaly, between -inf and inf, [radians] or [deg]
    """
    if degrees:
        nu_rad = np.deg2rad(nu)
        F_rad = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(nu_rad / 2))
        F = np.rad2deg(F_rad)
    else:    
        F = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(nu / 2))
    return F   

def E_to_Me(E,ecc,degrees=True):
    """
    Transform to Mean Anomaly from Eccentric Anomaly.

    Usage:
        >>> M = E_to_Me(E, ecc)
    Inputs:
        E -> [array-like,float] Eccentric Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        Me -> [array-like,float] Mean anomaly, [radians] or [deg]
    """
    twopi = Const.twopi
    if degrees:
        E_rad = np.deg2rad(E)
        Me_rad = E_rad - ecc * np.sin(E_rad)
        Me = np.rad2deg(Me_rad) % 360
    else:
        Me = (E - ecc * np.sin(E)) % twopi
    return Me  

def F_to_Mh(F,ecc,degrees=True):
    """
    Transform to Mean Anomaly from Hyperbolic Anomaly.

    Usage:
        >>> Mh = F_to_Mh(F, ecc)
    Inputs:
        F -> [array-like,float] Hyperbolic Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        Mh -> [array-like,float] Mean anomaly, [radians] or [deg]
    """
    twopi = Const.twopi
    if degrees:
        F_rad = np.deg2rad(F)
        Mh_rad = ecc * np.sinh(F_rad) - F_rad
        Mh = np.rad2deg(Mh_rad) % 360
    else:    
        Mh = (ecc * np.sinh(F) - F) % twopi
    return Mh

def nu_to_Me(nu,ecc,degrees=True):
    """
    Transform to Mean anomaly from True Anomaly for ellipse trajectories.

    Usage:
        >>> Me = nu_to_Me(nu, ecc)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        Me -> [array-like,float] Mean Anomaly, [radians] or [deg]
    """
    E = nu_to_E(nu,ecc,degrees)
    Me = E_to_Me(E,ecc,degrees)

    return Me  

def nu_to_Mp(nu,degrees=True):
    """
    Transform to Mean anomaly from True Anomaly for parabolic trajectories.

    Usage:
        >>> Mp = nu_to_Mp(nu)
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc = 1
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        Mp -> [array-like,float] Mean Anomaly, [radians] or [deg]
    """
    twopi = Const.twopi
    if degrees:
        nu_rad = np.deg2rad(nu)
        Mp_rad = np.tan(nu_rad/2)/2 + np.tan(nu_rad/2)**3/6
        Mp = np.rad2deg(Mp_rad) % 360
    else:    
        Mp = (np.tan(nu/2)/2 + np.tan(nu/2)**3/6) % twopi

    return Mp        

def nu_to_Mh(nu,ecc,degrees=True):
    """
    Transform to Mean anomaly from True Anomaly for hyperbolic trahectories.

    Usage:
        >>> Mh = nu_to_Mh(nu, ecc)
    
    Inputs:
        nu -> [array-like,float] True Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        Mh -> [array-like,float] Mean Anomaly, [radians] or [deg]
    """
    F = nu_to_F(nu,ecc,degrees)
    Mh = F_to_Mh(F,ecc,degrees)

    return Mh      

def Me_to_E(Me,ecc,degrees=True):
    """
    Transform to Eccentric Anomaly from Mean Anomaly.

    Usage:
        >>> E = Me_to_E(Me, ecc)
    Inputs:
        Me -> [array-like,float] Mean anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1 
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        E -> [array-like,float] Eccentric anomaly, [radians] or [deg]
    Notes:
        This uses a Newton iteration on the Kepler Equation.
    """
    if degrees: 
        Me = np.deg2rad(Me)
    else:
        Me = np.array(Me)  

    ecc = np.array(ecc)     

    if ecc.ndim ==0 and Me.ndim == 1:
        ecc = np.full_like(Me,ecc,dtype=float)
    if Me.ndim == 1:
        sols = []
        for M_i,ecc_i in zip(Me,ecc):
            sol = root_scalar(_kepler_equation_e,x0=M_i,args=(M_i,ecc_i), fprime=True, method='newton')
            sols.append(sol.root)
        E = np.array(sols)
    elif Me.ndim == 0:
        sol = root_scalar(_kepler_equation_e,x0=Me,args=(Me,ecc), fprime=True, method='newton')
        E = sol.root  
    if degrees: E = np.rad2deg(E)   
    return E  

def Mh_to_F(Mh,ecc,degrees=True):
    """
    Transform to Hyperbolic Anomaly from Mean Anomaly.

    Usage:
        >>> F = Mh_to_F(Mh, ecc)
    Inputs:
        Mh -> [array-like,float] Mean anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        F -> [array-like,float] Hyperbolic anomaly
    Notes:
        This uses a Newton iteration on the Kepler Equation.
    """
    if degrees: 
        Mh = np.deg2rad(Mh)
    else:
        Mh = np.array(Mh)  

    ecc = np.array(ecc)    
    if ecc.ndim ==0 and Mh.ndim == 1:
        ecc = np.full_like(Mh,ecc,dtype=float)  

    if Mh.ndim == 1:
        sols = []
        for M_i,ecc_i in zip(Mh,ecc):
            sol = root_scalar(_kepler_equation_h,x0=M_i,args=(M_i,ecc_i), fprime=True, method='newton')
            sols.append(sol.root)
        F = np.array(sols)
    elif Mh.ndim == 0:
        sol = root_scalar(_kepler_equation_h,x0=Mh,args=(Mh,ecc), fprime=True, method='newton')
        F = sol.root 
    if degrees: F = np.rad2deg(F)   
    return F     

def E_to_nu(E,ecc,degrees=True):
    """
    Transform to True anomaly from Eccentric Anomaly.

    Usage:
        >>> nu = E_to_nu(E, ecc)
    Inputs:
        E -> [array-like,float] Eccentric Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with 0 < ecc < 1 
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) or (-360,360)
    """
    twopi = Const.twopi
    if degrees:
        E_rad = np.deg2rad(E)
        nu_rad = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E_rad / 2))
        nu = np.rad2deg(nu_rad) % 360
    else:    
        nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2)) % twopi
    return nu 

def F_to_nu(F, ecc,degrees=True):
    """
    Transform to True anomaly from Hyperbolic Anomaly.

    Usage:
        >>> nu = F_to_nu(F, ecc)
    Inputs:
        F -> [array-like,float] Hyperbolic Anomaly, [radians] or [deg]
        ecc -> [array-like,float] Eccentricity with ecc > 1
        degrees -> [bool,optional,default=True] unit for angular variables 
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
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
    """
    E = Me_to_E(Me,ecc,degrees)
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
        degrees -> [bool,optional,default=True] unit for angular variables 
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
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
    """
    F = Mh_to_F(Mh,ecc,degrees)
    nu = F_to_nu(F,ecc,degrees) 
    return nu        

def coe2nse(a,ecc,inc,raan,argp,nu,degrees=True):
    """
    Convert to non-singular orbital elements from classical orbital elements for elliptic trajectories.
    The non-singular orbital elements exhibit no singularity for near-circular orbit, also known as the first kind of non-singular orbital elements. 

    Usage: 
        >>> a, inc, raan, xi, eta, l = coe2nse(a, ecc, inc, raan, argp, nu)
    Inputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        a -> [array-like,float] Semi-major axis
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        xi -> [array-like,float] non-singular parameter, ξ = e * cosω
        eta -> [array-like,float] non-singular parameter,η = e * cosω
        l -> [array-like,float] non-singular parameter, [radians] or [deg] , l = ω + M 
    """
    if ecc < 1:
        M = nu_to_Me(nu,ecc,degrees)
    else:
        raise Exception('Convertion to non-singular orbital elements only applies to elliptic trajectories.')    

    if degrees: 
        argp_rad = np.deg2rad(argp)
        xi = ecc * np.cos(argp_rad)
        eta = ecc * np.sin(argp_rad)
    else:
        xi = ecc * np.cos(argp)
        eta = ecc * np.sin(argp)
    l = argp + M

    nse = (a, inc, raan, xi, eta, l)

    return nse

def nse2coe(a,inc,raan,xi,eta,l,degrees=True):
    """
    Convert to classical orbital elements from non-singular orbital elements for elliptic trajectories.
    The non-singular orbital elements exhibit no singularity for near-circular orbit, also known as the first kind of non-singular orbital elements.

    Usage:
        >>> a, ecc, inc, raan, argp, nu = nse2coe(a, inc, raan, xi, eta, l)
    Inputs:
        a -> [array-like,float] Semi-major axis
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        xi -> [array-like,float] non-singular parameter, ξ = e * cosω
        eta -> [array-like,float] non-singular parameter, η = e * cosω
        l -> [array-like,float] non-singular parameter, [radians] or [deg], l = ω + M   
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
    """
    ecc = np.sqrt(xi ** 2 + eta ** 2)

    if not ecc < 1:
        raise Exception('Convertion from non-singular orbital elements only applies to elliptic trajectories.')   

    argp = np.arctan2(eta, xi)
    if degrees:
        argp = np.rad2deg(argp) % 360
        M = (l - argp)%360
    else:
        M = (l - argp) % Const.twopi

    nu = Me_to_nu(M,ecc,degrees) 
    coe = (a, ecc, inc,raan,argp, nu)   
        
    return coe

def coe2mee(a,ecc,inc,raan,argp,nu,degrees=True):
    """
    Convert to modified equinoctial orbital elements from classical orbital elements for elliptic trajectories. 
    The modified equinoctial orbital elements exhibit no singularity for near-circular orbit with inclinations close to 0 degrees. 
    It is also known as the second kind of non-singular orbital elements.

    Usage:
        >>> p, f, g, h, k, L = coe2mee(a, ecc, inc, raan, argp, nu)
    Inputs:
        a -> [array-like,float] Semi-major axis
        ecc -> [array-like,float] Eccentricity
        inc -> [array-like,float] Inclination, [radians] or [deg]
        omega -> [array-like,float] Longitude of ascending node, [radians] or [deg]
        argp -> [array-like,float] Argument of perigee, [radians] or [deg]
        nu -> [array-like,float] True anomaly, [radians] or [deg]
        degrees -> [bool,optional,default=True] unit for angular variables 
    Outputs:
        p -> [array-like,float] Semi-latus rectum, p = a*(1-e**2)
        f -> [array-like,float] Equinoctial parameter f, f = e * cos(Ω + ω)
        g -> [array-like,float] Equinoctial parameter g, g = e * sin(Ω + ω)
        h -> [array-like,float] Equinoctial parameter h, h = tan(i/2) * cos(Ω)
        k -> [array-like,float] Equinoctial parameter k, k = tan(i/2) * sin(Ω)
        L -> [array-like,float] Longitude, [radians] or [deg], L = Ω + ω + M
    """
    lonper = raan + argp
    p = a*(1-ecc**2)

    if ecc < 1:
        M = nu_to_Me(nu,ecc,degrees)
    else:
        raise Exception('Convertion to non-singular orbital elements only applies to elliptic trajectories.')  

    L = lonper + M

    if degrees:
        lonper = np.deg2rad(lonper)
        raan = np.deg2rad(raan)
        inc = np.deg2rad(inc)
    
    f = ecc * np.cos(lonper)
    g = ecc * np.sin(lonper)
    h = np.tan(inc / 2) * np.cos(raan)
    k = np.tan(inc / 2) * np.sin(raan)

    mee = (p, f, g, h, k, L)

    return mee

def mee2coe(p,f,g,h,k,L,degrees=True):
    """
    Convert to classical orbital elements from modified equinoctial orbital elements for elliptic trajectories.

    Usage:
        >>> a, ecc, inc, raan, argp, nu = mee2coe(p, f, g, h, k, L)
    Inputs:
        p -> [array-like,float] Semi-latus rectum, p = a*(1-e**2)
        f -> [array-like,float] Equinoctial parameter f, f = e * cos(Ω + ω)
        g -> [array-like,float] Equinoctial parameter g, g = e * sin(Ω + ω)
        h -> [array-like,float] Equinoctial parameter h, h = tan(i/2) * cos(Ω)
        k -> [array-like,float] Equinoctial parameter k, k = tan(i/2) * sin(Ω)
        L -> [array-like,float] Longitude, [radians] or [deg], L = Ω + ω + M
        degrees -> [bool,optional,default=True] unit for angular variables 
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

    if not ecc < 1:
        raise Exception('Convertion from non-singular orbital elements only applies to elliptic trajectories.') 

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
            
    nu = Me_to_nu(M,ecc,degrees)

    coe = (a, ecc, inc, raan, argp, nu)

    return coe 

def mean2osculating(mean_ele,epoch,meanref='TEME',oscuref='TEME',degrees=True):
    """
    Convert mean orbital elements to osculating orbital elements using sgp4/sdp4.

    Usage:
        >>> from astropy.time import Time
        >>> mean_ele = [7000,0.01,50,100,30,210] # in form of [a, e, i, Ω, ω, v]
        >>> epoch = Time('2022-06-07T08:09:12.345')
        >>> oscu_ele = mean2osculating(mean_ele,epoch)
    Inputs:
        mean_ele -> [list or array of float] mean elements for sgp4/sdp4
        epoch -> Object of class Astropy Time
        meanref -> [str,optional,default='TEME'] reference frame bound by the mean elements
        oscuref -> [str,optional,default='TEME'] reference frame bound by the osculating elements
        degrees -> [bool,optional,default=True] unit of the angular variable of orbital elements
    Outputs:
        oscu_ele -> [list or array of float] osculating elements for sgp4/sdp4
    """
    ts = data_prepare.ts
    t = ts.from_astropy(epoch)
    
    if meanref in ['ECI','J2000','ICRF']:
        gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(epoch)
        mean_ele = coe_trans(gcrf2teme_mat,mean_ele,degrees)   
        
    a = mean_ele[0]/Const.Re_sgp4
    ecc = mean_ele[1]

    if degrees: inc,raan,argp,M = np.deg2rad(mean_ele[2:])

    satrec = sgp4init(a,ecc,inc,raan,argp,M,epoch)
    sat = EarthSatellite.from_satrec(satrec, ts)   
        
    if oscuref == 'TEME':
        sat_rv = np.hstack(sat._position_and_velocity_TEME_km(t)[:2]) 
    elif oscuref in ['ECI','J2000','ICRF']:
        sat_t = sat.at(t)
        sat_rv = np.hstack([sat_t.xyz.km,sat_t.velocity.km_per_s]) 

    osculating_nu = rv2coe(sat_rv,Const.mu_sgp4,degrees) # WGS72
    osculating_Me = osculating_nu.copy()
    osculating_Me[-1] = nu_to_Me(osculating_nu[-1],osculating_nu[1],degrees)  
    
    return osculating_Me

def func_oscu2mean(mean_ele,oscu_ele,epoch,degrees):
    """
    Auxiliary functions converting osculating orbital elements to mean orbital elements using sgp4/sdp4.

    Inputs:
        mean_ele -> [list or array of float] mean elements for sgp4/sdp4
        oscu_ele -> [list or array of float] osculating elements for sgp4/sdp4
        epoch -> Object of class Astropy Time
        degrees -> [bool] unit of the angular variable of orbital elements
    Outputs:
        deviation -> [float] O-C of state vector
    """
    ts = data_prepare.ts
    ta = Time(epoch)
    t = ts.from_astropy(ta)

    a,ecc,inc,raan,argp,M = mean_ele 

    satrec = sgp4init(a,ecc,inc,raan,argp,M,epoch)
    sat = EarthSatellite.from_satrec(satrec, ts)

    sat_rv = np.hstack(sat._position_and_velocity_TEME_km(t)[:2])  
    sat_rv_nd = np.hstack([sat_rv[:3]/Const.Re_sgp4,sat_rv[3:]/Const.v_sgp4])

    oscu_ele_nu = oscu_ele.copy()

    oscu_ele_nu[-1] = Me_to_nu(oscu_ele[-1],oscu_ele[1],degrees)
    oscu_rv = coe2rv(oscu_ele_nu,Const.mu_sgp4,degrees)
    oscu_rv_nd = np.hstack([oscu_rv[:3]/Const.Re_sgp4,oscu_rv[3:]/Const.v_sgp4])

    return norm(sat_rv_nd - oscu_rv_nd)    

def osculating2mean(oscu_ele,epoch,oscuref='TEME',meanref='TEME',degrees=True):
    """
    Convert osculating orbital elements to mean orbital elements using sgp4/sdp4.

    Usage:
        >>> from astropy.time import Time
        >>> oscu_ele = [7000,0.01,50,100,30,210] # in form of [a, e, i, Ω, ω, v]
        >>> epoch = Time('2022-06-07T08:09:12.345')
        >>> mean_ele = osculating2mean(oscu_ele,epoch)
    Inputs:
        oscu_ele -> [list or array of float] osculating elements for sgp4/sdp4
        epoch -> Object of class Astropy Time
        oscuref -> [str,optional,default='TEME'] reference frame bound by the osculating elements
        meanref -> [str,optional,default='TEME'] reference frame bound by the mean elements
        degrees -> [bool,optional,default=True] unit of the angular variable of orbital elements
    Outputs:
        mean_ele -> [list or array of float] mean elements for sgp4/sdp4
    """
    twopi = Const.twopi

    if oscuref in ['ECI','J2000','ICRF']:
        gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(Time(epoch))
        oscu_ele = coe_trans(gcrf2teme_mat,oscu_ele,degrees)  

    oscu_ele_a = oscu_ele[0]/Const.Re_sgp4
    oscu_ele_e = oscu_ele[1]

    if degrees: 
        oscu_ele_ = np.deg2rad(oscu_ele[2:])
    else:
        oscu_ele_ = oscu_ele[2:]
            
    mean_ele = np.hstack([oscu_ele_a,oscu_ele_e,oscu_ele_])

    bounds = [(oscu_ele_a*0.9,oscu_ele_a*1.1),(oscu_ele_e*0.8,oscu_ele_e*1.2),(0,np.pi),(0,twopi),(0,twopi),(0,twopi)]

    res = minimize(func_oscu2mean,mean_ele,args = (oscu_ele,epoch,degrees),bounds=bounds)
    if not res.success: raise Exception(res.message)
    mean_ele = res.x

    if degrees: 
        mean_ele_ = np.rad2deg(mean_ele[2:])
    else:
        mean_ele_ = mean_ele[2:]

    mean_ele = np.hstack([mean_ele[0]*Const.Re_sgp4,mean_ele[1],mean_ele_])

    if meanref in ['ECI','J2000','ICRF']:
        gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(Time(epoch))
        mean_ele = coe_trans(teme2gcrf_mat,mean_ele,degrees)  

    return mean_ele    

def coe_trans(trans_matrix,coe_from,degrees=True):
    """
    Converting classical orbital elements between two reference frames, especially 'TEME' and 'ICRF'. The reference frame must be defined as hand-right.

    Usage:
        >>> coe_to = coe_trans(trans_matrix,coe_from)
    Inputs:
        trans_matrix -> [nx3x3 array-like] single or mutiple transformation matrix
        coe_from -> [nx6] classical orbital elements in source reference frame
        degrees -> [bool,optional,default=True] unit of the angular variable of orbital elements
    Outputs:
        coe_to -> [nx6] classical orbital elements in target reference frame  
    """
    twopi = Const.twopi
    a,ecc,inc,raan,argp,M = np.array(coe_from).T
    X,Y,Z = euler2vectors(np.array([raan,inc,argp]).T,degrees)

    hx_,hy_,hz_ = h_uec_ = Matrix_dot_Vector(trans_matrix,Z).T 
    x_,y_,z_ = xyz_ = Matrix_dot_Vector(trans_matrix,X).T
        
    inc_ = np.arccos(np.clip(hz_,-1,1))

    # for inc_ != [0,pi]
    raan_ = np.arctan2(hx_,-hy_)
    cos_argp_ = x_*np.cos(raan_) + y_*np.sin(raan_)
    sin_argp_ = (y_*np.cos(raan_) - x_*np.sin(raan_))/np.cos(inc_)
    argp_ = np.arctan2(sin_argp_,cos_argp_)

    if np.isscalar(inc):
        if inc in [0,np.pi]:
            raan_ = 0
            argp_ = np.arctan2(y_,x_)
    else:    
        # for inc_ == [0,pi]
        inc_flag = np.in1d(inc_,[0,np.pi])
        raan_[inc_flag] = 0
        argp_[inc_flag] = np.arctan2(y_,x_)[inc_flag]

    raan_,argp_ = raan_%twopi,argp_%twopi
    if degrees: inc_,raan_,argp_ = np.rad2deg([inc_,raan_,argp_])
    ele_ = np.array([a,ecc,inc_,raan_,argp_,M]).T

    return ele_    

def rv_trans(trans_matrix,rv_from):
    """
    Converting orbital state vector between two reference frames.

    Usage:
        >>> rv_to = rv_trans(trans_matrix,rv_from)
    Inputs:
        trans_matrix -> [nx3x3 array-like] single or mutiple transformation matrix
        rv_from -> [nx6] orbital state vector in source reference frame
    Outputs:
        rv_to -> [nx6] orbital state vector in target reference frame  
    """
    rv_from_T = np.array(rv_from).T
    xyz_from,vxyz_from = rv_from_T[:3].T,rv_from_T[3:].T
    xyz_to = Matrix_dot_Vector(trans_matrix,xyz_from)   
    vxyz_to = Matrix_dot_Vector(trans_matrix,vxyz_from)
    rv_to = np.hstack([xyz_to,vxyz_to])
    return rv_to          