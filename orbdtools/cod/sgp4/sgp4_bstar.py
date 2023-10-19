import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar
from skyfield.api import EarthSatellite

from .sgp4_init import no_kozai_calculate,sgp4init_no_kozai
from ...transform.orbele_trans import Me_to_nu,coe_trans
from ...transform.kep_rv_trans import coe2rv
from ...transform.frame_trans import gcrf_teme_mat
from ...utils import Const,data_prepare

def func_bstar(bstar,no_kozai1,ecc1,inc1,raan1,argp1,M1,epoch1,ele2_o_rv,epoch2,ts):
    """
    Residual function for estimation of bstar.
    """
            
    satrec1 = sgp4init_no_kozai(no_kozai1,ecc1,inc1,raan1,argp1,M1,epoch1,bstar=bstar)
    sat = EarthSatellite.from_satrec(satrec1,ts)
    
    t2 = ts.from_astropy(epoch2)
    sat_t = sat.at(t2)
    satm = sat.model
    ele2_teme_c = np.array([satm.am,satm.em,satm.im,satm.Om,satm.om,satm.mm])
    ele2_c_rv = coe2rv(ele2_teme_c,1,False)
    
    return norm(ele2_o_rv - ele2_c_rv)  

def bstar_estimate(ele1,ele2,epoch1,epoch2,degrees=True,ref='TEME'):
    """
    Estimate the Bstar involved in SGP4 propagator based on Mean Elements at two epoch times.

    Usage:
        >>> from astropy.time import Time
        >>> from orbdtools.cod.sgp4 import sgp4_bstar
        >>> # Mean Elements in TEME at epoch1 in form of [a, e, i, Ω, ω, M]
        >>> epoch1 = Time('2022-05-23T22:00:02.000Z')
        >>> ele1_teme = [1.071459,0.000154,53.2175,182.0098,49.6208,205.5672]
        >>> # Mean Elements in TEME at epoch2 in form of [a, e, i, Ω, ω, M]
        >>> epoch2 = Time('2022-05-25T08:38:34.000Z')
        >>> ele2_teme = [1.072465,0.000159,53.2175,175.255990,58.766942,261.018083]
        >>> # Estimate B*
        >>> bstar = sgp4_bstar.bstar_estimate(a_ecc_inc1,a_ecc_inc2,epoch1,epoch2,degrees=True)
        >>> print('Estimated B*: ',bstar)
        >>> # Estimated B*:  -0.02128338135900665
    Inputs:
        ele1 -> [array of float] Mean Elements at the first epoch in form of [a,ecc,inc,raan,argp,M], where the unit of 'a' is the equatorial radius of the earth defined by WGS72(6378.135km).
        ele2 -> [array of float] Mean Elements at the second epoch in form of [a,ecc,inc,raan,argp,M]
        epoch1 -> [Astropy Time] The first epoch
        epoch2 -> [Astropy Time] The second epoch
        degrees -> [bool,optional,default=True] Unit of angular variables in orbital elements. If True, angular variables are in [deg], otherwise in [rad].
        ref -> [str,optional,default='TEME'] Reference frame in which the Mean Elements defined. Available options include 'TEME', and 'GCRF'/'J2000'/'ECI'. 
    Outputs:
        bstar -> [float] Value of bstar in unit of 1/radius
    """
    ts = data_prepare.ts

    if ref == 'TEME':
        ele1_teme,ele2_teme = ele1,ele2
    elif ref in ['GCRF','J2000','ECI']:
        # convert to teme
        gcrf2teme_mat1,teme2gcrf_mat1 = gcrf_teme_mat(epoch1)
        gcrf2teme_mat2,teme2gcrf_mat2 = gcrf_teme_mat(epoch2)
        ele1_teme = coe_trans(gcrf2teme_mat1,ele1,degrees)
        ele2_teme = coe_trans(gcrf2teme_mat2,ele2,degrees)
    else:
        raise Exception("Unsupported reference frame. Available options include 'TEME' and 'GCRF'/'J2000'/'ECI'.")    

    a1,ecc1,inc1,raan1,argp1,M1 = ele1_teme 
    a2,ecc2,inc2,raan2,argp2,M2 = ele2_teme   

    nu2 = Me_to_nu(M2,ecc2,degrees)
    ele2_o_rv = coe2rv([a2,ecc2,inc2,raan2,argp2,nu2],1,degrees) 

    if degrees: inc1,raan1,argp1,M1 = np.deg2rad([inc1,raan1,argp1,M1])
    no_kozai1 = no_kozai_calculate(a1,ecc1,inc1)
    
    res = minimize_scalar(func_bstar,args = (no_kozai1,ecc1,inc1,raan1,argp1,M1,epoch1,ele2_o_rv,epoch2,ts))
    if not res.success: raise Exception(res.message)
    bstar = res.x
    return bstar   

def sgp4_C2(no_kozai,ecc,inc):
    """
    Calculate C2 defined in the SGP4 model.

    Usage:
        >>> C2,a0dp = sgp4_C2(no_kozai,ecc,inc)
    Inputs:
        no_kozai -> [float] Mean motion with unit of [rad/min] in Kozai theory
        ecc -> [float] Eccentricity
        inc -> [float] Inclination in [rad]
    Outputs:
        C2 -> [float] C2 defined in the SGP4 model with unit of [Length unit/Time unit], 
        where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km, 
        and [Time unit] is equal to ([Length unit]^3/[mu])**0.5 with mu = 398600.8 [km^3/s^2] 
        a0dp -> [float] Original semimajor with unit of [Length unit] in Brouwer theory, as opposed to that in Kozai theory.
    Reference:
        Hoots F R. Spacetrack report no. 3, models for propagation of norad element sets[J]. http://www. itc. nl/-bakker/orbit. html, 1980.    
    """
    ke = Const.ke_sgp4
    k2 = Const.k2_sgp4
    Re = Const.Re_sgp4
    
    no_kozai = no_kozai/60*Const.T_sgp4 # from [rad/min] to [rad/unit time]
    
    a1 = (ke / no_kozai)**(2/3)
    
    theta = cosi0 = np.cos(inc)
    x3thm1 = 3 * cosi0**2 - 1
    
    del1 = 1.5 * k2 * x3thm1 / a1**2 / (1 - ecc**2)**1.5
    a0 = a1*(1 - del1/3 - del1**2 - 134/81*del1**3)
    del0 = 1.5 * k2 * x3thm1 / a0**2 / (1 - ecc**2)**1.5
    
    n0dp = no_kozai / (1 + del0)
    a0dp = a0 / (1 - del0)
    
    s =  1 + 78 / Re
    q_sub_s = (120 - 78) / Re
    
    perige = (a0dp * (1 - ecc) - 1) * Re # in km
    
    if perige < 98:
        s_star = 1 + 20/Re
    elif 98 <= perige <= 156:
        s_star = a0dp*(1 - ecc) - s + 1
    else:
        s_star = s
        
    q_sub_s_star = (q_sub_s + s - s_star)**4
    
    tsi = 1 / (a0dp - s_star)
    
    beta0 = np.sqrt(1 - ecc**2)
    eta = a0dp*ecc*tsi
    
    C2 = q_sub_s_star*tsi**4*n0dp*(1-eta**2)**(-3.5)*(a0dp*(1+1.5*eta**2+4*ecc*eta+ecc*eta**3)+1.5*k2*tsi/(1 - eta**2)*(-0.5+1.5*theta**2)*(8+24*eta**2+3*eta**4))
    
    return C2,a0dp

def bstar_calculate(a_ecc_inc1,a_ecc_inc2,epoch1,epoch2,degrees=True):
    """
    Calculate the Bstar involved in SGP4 propagator based on Mean Elements at two epoch times.

    Usage:
        >>> from astropy.time import Time
        >>> from orbdtools.cod.sgp4 import sgp4_bstar
        >>> # Mean Elements in TEME at epoch1 in form of [a, e, i]
        >>> epoch1 = Time('2022-05-23T22:00:02.000Z')
        >>> a_ecc_inc1 = [1.071459,0.000154,53.2175]
        >>> # Mean Elements in TEME at epoch2 in form of [a, e, i]
        >>> epoch2 = Time('2022-05-25T08:38:34.000Z')
        >>> a_ecc_inc2 = [1.072465,0.000159,53.2175]
        >>> # Calculate B*
        >>> bstar = sgp4_bstar.bstar_calculate(a_ecc_inc1,a_ecc_inc2,epoch1,epoch2,degrees=True)
        >>> print('Calculated B*: ',bstar)
        >>> # Calculated B*:  -0.02118445232232004
    Inputs:
        a_ecc_inc1 -> [array of float] Mean Elements in TEME at the first epoch in form of [a,ecc,inc], where the unit of 'a' is the equatorial radius of the earth defined by WGS72(6378.135km).
        a_ecc_inc2 -> [array of float] Mean Elements in TEME at the second epoch in form of [a,ecc,inc].
        epoch1 -> [Astropy Time] The first epoch
        epoch2 -> [Astropy Time] The second epoch
        degrees -> [bool,optional,default=True] Unit of angular variables in orbital elements. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        bstar -> [float] Value of bstar in unit of 1/radius
    """
    inc1 = a_ecc_inc1[2]
    inc2 = a_ecc_inc2[2]  

    if degrees: inc1,inc2 = np.deg2rad([inc1,inc2]) 

    a1,ecc1 = a_ecc_inc1[:2]
    no_kozai1 = no_kozai_calculate(a1,ecc1,inc1)
    C2_1,_a0dp1 = sgp4_C2(no_kozai1,ecc1,inc1)

    a2,ecc2 = a_ecc_inc2[:2]
    no_kozai2 = no_kozai_calculate(a2,ecc2,inc2)
    C2_2,_a0dp2 = sgp4_C2(no_kozai2,ecc2,inc2)

    C2 = np.average([C2_1,C2_2])
    n1,n2 = np.sqrt(1/a1**3),np.sqrt(1/a2**3)
    tof = (epoch2 - epoch1).sec/Const.T_sgp4
    XNDT20 = (n2 - n1)/2/tof
    n = np.average([n1,n2])
    bstar = XNDT20/n/C2/1.5
    
    return bstar