import numpy as np
from astropy.coordinates import spherical_to_cartesian
from astropy import units
from loess.loess_1d import loess_1d

from . import Const  

def get_mid_point(t,*args):
    """
    Extract the median point from the observation data.

    Usage:
        >>> t_mid,xyz_site_mid,radec_mid = get_mid_point(t,xyz_site,radec)
    Inputs:
        t -> [array-like,Astropy Time or Skyfield Time] time
        xyz_site -> [array like, optional] cartesian coordinates of site
        radec -> [array like, optional] RA and Dec of space object
    Outputs:
        t_mid -> [object] median time in form of Astropy Time or Skyfield Time
        xyz_site_mid -> [array like, optional] cartesian coordinates corresponding to the median time
        radec_mid -> [array like, optional] RA and Dec corresponding to the median time  
    """
    tof = t[-1] - t[0]
    indices = np.zeros(3,dtype=int)

    mt = t[0] + tof/2 # node moment of arc
    mi = np.argmin(np.abs(t - mt)) # index corresponding to node moment

    out = [t[mi]]       

    for cont in args:
        out.append(cont[mi])    

    return out   

def lowess_smooth_optical(ta,radec,frac=0.5,T_scale=Const.T_nd_unit):
    """
    Remove outliers in optical data with the method of LOWESS (Locally Weighted Scatterplot Smoothing)
    Here, LOWESS uses a weighted **linear regression** in default.

    Usage:
        >>> flag = lowess_smooth_optical(ta,radec)
    Inputs:
        ta -> [array of Astropy Time] Time sequence
        radec -> [2D array] RA and Dec of space object, [deg]
        frac -> [float,optional,default=0.5] The fraction of the data used in local regression. The value of fraction is between 0 and 1.
        T_scale -> [float,optional,default=Const.T_nd_unit] Time scale, defaults to sqrt([L_nd]**3/[mu_nd]), where [L_nd]=6378.137 [km] and [mu_nd]=398600.4418 [km^3/s^2] are the equatorial radius and GM of the Reference Earth Model - WGS84 respectively.
    Outputs:
        flag -> [array of bool] Flag of data points. If False, the data point is an outlier.        
    """
    # Calculate line of sight
    t_ = (ta - ta[0]).sec/T_scale
    los = spherical_to_cartesian(1,radec[:,1]*units.deg,radec[:,0]*units.deg) 
    los = np.stack([los_i.value for los_i in los]).T  

    _, los_x, w_los_x = loess_1d(t_,los[:,0],frac=frac)
    _, los_y, w_los_y = loess_1d(t_,los[:,1],frac=frac)
    _, los_z, w_los_z = loess_1d(t_,los[:,2],frac=frac)

    w = 0.34 # This can be used to identify outliers: wout < 0.34 for outliers with deviations > 4sigma.

    flag = (w_los_x > w) & (w_los_y > w) & (w_los_z > w) 

    return flag

def lowess_smooth_radar(ta,azalt,r,frac=0.5,L_scale=Const.L_nd_unit,T_scale=Const.T_nd_unit):
    """
    Remove outliers in radar data with the method of LOWESS (Locally Weighted Scatterplot Smoothing)
    Here, LOWESS uses a weighted **linear regression** in default.

    Usage:
        >>> flag = lowess_smooth_radar(ta,radec,r)
    Inputs:
        ta -> [array of Astropy Time] Time sequence
        azalt -> [2D array] Az and Alt of space object, [deg]
        r -> [array] Slant distance of the space object relative to the site, [km]
        frac -> [float,optional,default=0.5] The fraction of the data used in local regression. The value of fraction is between 0 and 1.  
        L_scale -> [float,optional,default=Const.L_nd_unit] Length scale, defaults to the equatorial radius of the Reference Earth Model - WGS84.  
        T_scale -> [float,optional,default=Const.T_nd_unit] Time scale, defaults to sqrt([L_nd]**3/[mu_nd]), where [L_nd]=6378.137 [km] and [mu_nd]=398600.4418 [km^3/s^2] are the equatorial radius and GM of the Reference Earth Model - WGS84 respectively.
    Outputs:
        flag -> [array of bool] Flag of data points. If False, the data point is an outlier.        
    """
    # Calculate the cartesian coordinates of space object relative to site
    t_ = (ta - ta[0]).sec/T_scale
    r_nd = r/L_scale
    xyz = spherical_to_cartesian(r_nd,azalt[:,1]*units.deg,azalt[:,0]*units.deg) 
    xyz = np.stack([xyz_i.value for xyz_i in xyz]).T  

    _, x, w_x = loess_1d(t_,xyz[:,0],frac=frac)
    _, y, w_y = loess_1d(t_,xyz[:,1],frac=frac)
    _, z, w_z = loess_1d(t_,xyz[:,2],frac=frac)

    w = 0.34 # This can be used to identify outliers: wout < 0.34 for outliers with deviations > 4sigma.

    flag = (w_x > w) & (w_y > w) & (w_z > w) 

    return flag   