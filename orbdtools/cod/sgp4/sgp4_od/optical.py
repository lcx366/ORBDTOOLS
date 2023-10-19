import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares
from skyfield.api import EarthSatellite

from ..sgp4_init import sgp4init_no_kozai,export_tle,no_kozai_calculate
from ..sgp4_bstar import sgp4_C2
from ....transform.frame_trans import gcrf_teme_mat
from ....transform.orbele_trans import coe_trans,Me_to_nu,nu_to_Me
from ....utils import Const,data_prepare

def fun_resi_optical(ele_bstar,ta0,t,R_vec,losnp):
    """
    Residual function for least squares estimation of orbit elements and bstar using optical angle measurements with SGP4 propagator.

    Usage:
        >>> residuals = fun_resi_optical(ele_bstar,ta0,t,R_vec,losnp)
    Inputs:
        ele_bstar -> [list of float] Orbit elements in form of [no_kozai,ecc,inc,raan,argp,M,bstar], where
            no_kozai -> [float] Mean motion in unit of [rad/min]
            ecc -> [float] Eccentricity
            inc -> [float] Inclination in [rad]
            raan -> [float] Right Ascension of Ascending Node in [rad]
            argp -> [float] Argument of Perigee in [rad]
            M -> [float] Mean Anomaly in [rad]
        ta0 -> [Object of Astropy Time] Epoch of the orbital elements in UTC
        t -> [Object of Skyfield Time] Time sequence from observations
        R_vec -> [2D array] Cartesian coordinates of site in GCRF(Geocentric Celestial Reference Frame) in in [Length unit], 
        where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km.
        losnp -> [2D array] Line of Sight(LOS) from site to space object in GCRF
    Outputs:
        residuals -> [Array of float] Difference between one and cosine of the angle spaned by the observational LOS and the calculated LOS at the Median epoch
    """
    ts = data_prepare.ts

    no_kozai,ecc,inc,raan,argp,M,bstar = ele_bstar

    satrec = sgp4init_no_kozai(no_kozai,ecc,inc,raan,argp,M,ta0,bstar=bstar) 
    sat = EarthSatellite.from_satrec(satrec, ts)
    sat_t = sat.at(t)

    rs_vec = sat_t.xyz.km.T/Const.Re_sgp4
    rhos_vec = rs_vec - R_vec

    proj_insight = np.sum(rhos_vec*losnp,axis=1)
    residuals = proj_insight/norm(rhos_vec,axis=1) - 1

    return residuals   

def sgp4_od_optical(ele0,ta0,ta,xyz_site,losnp,params):
    """
    Cataloging orbit determination using optical angle measurements with SGP4 propagator. 

    Usage:    
        >>> ele_dict,tle_str,rms = sgp4_od_optical(ele0,ta0,ta,xyz_site,losnp,params)
    Inputs:
        ele0 -> [list of float] Initial orbit elements in form of [a, e, i, Ω, ω, M], where
            a -> [float] Semimajor in [Length unit], where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km. 
            Considering that the orbit elements are only initial values, [Length unit] can take slightly different values, such as the equatorial radius defined by the WGS84 model.
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [deg]
            raan -> [float] Longitude of ascending node, [deg]
            argp -> [float] Argument of perigee, [deg]
            M -> [float] Mean anomaly, [deg]
        ta0 -> [Object of Astropy Time] Epoch of the orbital elements in UTC
        ta -> [Object of Astropy Time] Time sequence from observations
        xyz_site -> [2D array with shape of nx3] Cartesian coordinates of the site in GCRF(Geocentric Celestial Reference Frame), [km]
        losnp -> [2D array with shape of nx3] Line of Sight(LOS) from site to space object in GCRF
        params -> [list] Extended parameters for generating TLE:
            satnum -> [int] Satellite Catalog number, such as NORAD ID.
            reff -> [str] Reference Frame at which the initial orbit elements is defined. Available options include 'GCRF','J2000', 'ECI', and 'TEME'.
            Here, 'GCRF', 'J2000', and 'ECI' are treated as equivalent.
            bstar -> [float] The drag term, or radiation pressure coefficient in unit of [1/earth radii]. Its absolute value is usually less than 1.
            nddot -> [float] One sixth the second derivative of mean motion in unit of [rad/min^3]. Note that this item is only used for SGP, not SGP4. 
            classification -> [str] Classification (U: unclassified, C: classified, S: secret)
            intldesg -> [str] International Designator, where YY is the last two digits of launch year, XXX is the launch number of the year, and A is the piece of the launch.
            elnum -> [int] Element set number. Incremented when a new TLE is generated for this object.
            revnum -> [int] Revolution number at epoch 
    Outputs:
        ele_dict -> [dict] Improved orbit elements with estimated bstar and convergence state of the solution
        tle_str -> [tuple] Two line elements
        rms -> [float] RMS of O-C
    """
    twopi = Const.twopi

    ts = data_prepare.ts
    t = ts.from_astropy(ta)

    R_vec = xyz_site / Const.Re_sgp4 # Convert length to normalized unit

    satnum,reff,bstar,nddot,classification,intldesg,elnum,revnum = params

    if reff in ['GCRF','J2000','ECI']:
        # Convert J2000 to TEME
        gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta0)
        ele0 = coe_trans(gcrf2teme_mat,ele0)
    elif reff != 'TEME':
        raise Exception("Unsupported reference frame. Available options include 'TEME' and 'GCRF'/'J2000'/'ECI'.")    

    a,ecc = ele0[:2]
    inc,raan,argp,M = np.deg2rad(ele0[2:])
    no_kozai = no_kozai_calculate(a,ecc,inc)

    ele_bstar = [no_kozai,ecc,inc,raan,argp,M,bstar]
    res = least_squares(fun_resi_optical,ele_bstar, args=(ta0,t,R_vec,losnp),loss='huber',bounds=([0,0,0,0,0,0,-1],[1,1,np.pi,twopi,twopi,twopi,1]))
    no_kozai_,ecc_,inc_,raan_,argp_,M_,bstar_ = res.x

    O_C = res.fun
    rms = np.sqrt(np.dot(O_C,O_C)/len(O_C))

    nu_ = Me_to_nu(M_,ecc_,degrees=False)
    _C2,a_ = sgp4_C2(no_kozai_,ecc_,inc_)
    h_ = np.sqrt(a_ * (1 - ecc_**2))

    if res.success:
        status = 'success'
    else:
        status = 'failed' 

    n_ = np.sqrt(1/a_**3)
    # Calculate approximatelly the ballistic coefficient or one half the first time derivative of the mean motion in unit of [rad/min^2]. Note that this item is only used for SGP, not SGP4. 
    ndot_ = 1.5*bstar_*_C2*n_/(Const.T_sgp4/60)   
    
    tle_str = export_tle(no_kozai_,ecc_,inc_,raan_,argp_,M_,ta0,bstar_,ndot_,nddot,satnum,classification,intldesg,elnum,revnum)
  
    inc_,raan_,argp_,nu_,M_ = np.rad2deg(inc_),np.rad2deg(raan_),np.rad2deg(argp_),np.rad2deg(nu_),np.rad2deg(M_) 
    ele_dict = {'epoch':ta0.isot,'a':a_,'ecc':ecc_,'inc':inc_,'raan':raan_,'argp':argp_,'nu':nu_,'M':M_,'h':h_,'bstar':bstar_,'status':status}

    return ele_dict,tle_str,rms