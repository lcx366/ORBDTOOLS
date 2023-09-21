import numpy as np
from numpy.linalg import norm
from scipy.optimize import least_squares
from astropy.time import Time
from astropy.coordinates import spherical_to_cartesian
from astropy import units
from skyfield.api import EarthSatellite

from .sgp4_init import sgp4init_no_kozai,export_tle,no_kozai_calculate
from .sgp4_bstar import sgp4_C2
from ...transform.frame_trans import gcrf_teme_mat
from ...transform.orbele_trans import coe_trans,Me_to_nu,nu_to_Me
from ...utils import Const,data_prepare

def from_ele_dict(ele_dict):
    ele = [ele_dict['a'],ele_dict['ecc'],ele_dict['inc'],ele_dict['raan'],ele_dict['argp'],ele_dict['M']]
    t0 = ele_dict['epoch']

    return ele,t0

def fun_resi_optical(ele_bstar,ta0,t,R_vec,losnp,nd=True):
    '''
    residual, used for lsq
    '''

    ts = data_prepare.ts

    a,ecc,inc,raan,argp,M,bstar = ele_bstar

    satrec = sgp4init_no_kozai(a,ecc,inc,raan,argp,M,ta0,bstar=bstar) 
    sat = EarthSatellite.from_satrec(satrec, ts)
    sat_t = sat.at(t)

    rs_vec = sat_t.xyz.km.T/Const.Re_sgp4
    rhos_vec = rs_vec - R_vec
    proj_insight = np.sum(rhos_vec*losnp,axis=1)[:,None]*losnp
    residuals = norm(rhos_vec - proj_insight,axis=1)/norm(proj_insight,axis=1)

    return residuals  

def sgp4_od_optical(ele_ref,t,pos_obs,radec,**kwargs):
    """
    for ele_ref, degrees = True
    """

    twopi = Const.twopi

    ts = data_prepare.ts
    ta = Time(t)
    t = ts.from_astropy(ta)

    ele,t0 = from_ele_dict(ele_ref)

    ele_ref_keys = ele_ref.keys()

    if 'satnum' in kwargs: 
        satnum = kwargs['satnum'] 
    else:
        satnum = ele_ref['satnum']    

    if 'bstar' in kwargs: 
        bstar = kwargs['bstar']
    elif 'bstar' in ele_ref_keys:  
        bstar = ele_ref['bstar']
    else:    
        bstar = 0  

    if 'ndot' in kwargs: 
        ndot = kwargs['ndot']
    else:  
        ndot = 0  

    if 'nddot' in kwargs: 
        nddot = kwargs['nddot']
    else:  
        nddot = 0

    if 'classification' in kwargs: 
        classification = kwargs['classification']
    else:    
        classification = 'U'

    if 'intldesg' in kwargs: 
        intldesg = kwargs['intldesg']
    else:    
        intldesg = 'YYXXXA'

    if 'elnum' in kwargs: 
        elnum = kwargs['elnum']
    else:    
        elnum = 1

    if 'revnum' in kwargs: 
        revnum = kwargs['revnum']
    else:    
        revnum = 0

    # Convert J2000 to TEME
    ta0 = Time(t0)
    gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta0)
    ele = coe_trans(gcrf2teme_mat,ele)

    a,ecc = ele[:2]
    inc,raan,argp,M = np.deg2rad(ele[2:])
    no_kozai = no_kozai_calculate(a,ecc,inc)
    ele_bstar = [no_kozai,ecc,inc,raan,argp,M,bstar]

    R_vec = pos_obs / Const.Re_sgp4 # Convert length to normalized unit
    losnp = spherical_to_cartesian(1,radec[:,1]*units.deg,radec[:,0]*units.deg)
    losnp = np.stack([los_i.value for los_i in losnp]).T

    res = least_squares(fun_resi_optical,ele_bstar, args=(ta0,t,R_vec,losnp),loss='huber',bounds=([0,0,0,0,0,0,-1],[1,1,np.pi,twopi,twopi,twopi,1]))
    no_kozai_,ecc_,inc_,raan_,argp_,M_,bstar_ = res.x

    nu_ = Me_to_nu(M_,ecc_,degrees=False)
    _C2,a_ = sgp4_C2(no_kozai_,ecc_,inc_)
    h_ = np.sqrt(a_ * (1 - ecc_**2))

    if res.success:
        status = 'success'
    else:
        status = 'failed' 

    # calculate ndot approximatelly
    n_ = np.sqrt(1/a_**3)
    ndot_ = 1.5*bstar_*_C2*n_/(Const.T_sgp4/60)
    
    tle_string = export_tle(no_kozai_,ecc_,inc_,raan_,argp_,M_,ta0,bstar_,ndot_,nddot,satnum,classification,intldesg,elnum,revnum)
  
    inc_,raan_,argp_,nu_,M_ = np.rad2deg(inc_),np.rad2deg(raan_),np.rad2deg(argp_),np.rad2deg(nu_),np.rad2deg(M_) 

    ele_lsq_dict = {'epoch':t0,'a':a_,'ecc':ecc_,'inc':inc_,'raan':raan_,'argp':argp_,'nu':nu_,'M':M_,'h':h_,'bstar':bstar_,'status':status}

    return [ele_lsq_dict] ,tle_string 