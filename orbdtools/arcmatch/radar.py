import numpy as np
from numpy.linalg import norm,lstsq
from sgp4.api import SatrecArray
from astropy.time import Time
from astropy.coordinates import SkyCoord,spherical_to_cartesian,cartesian_to_spherical

from .parse_tle import load_tle_file
from ..transform.frame_trans import gcrf_teme_mat,ECI_RADAR_mat
from ..utils.preprocessing import get_mid_point
from ..utils.math import Matrix_dot_Vector

def arcsat_match(tle,ta,xyz_site,orbele_site,azalt,r,threshold_pre=[10,10],threshold_deep=[2000,5],threshold_slope=[400,2]):
    """
    Match the observation arc based on radra measurement data(range+angle) to space objects in TLE file.

    Usage:
        >>> code_match,satnum,disp_match = arcsat_match(tle,ta,xyz_site,radec,r)
    Inputs:
        tle -> object of class TLE
        ta -> [Astropy Time] time
        xyz_site -> [2D array] cartesian coordinates of site, [km]
        azalt -> [2D array] Az and Alt of space object, [deg]
        r -> [array] Slant distance of the space object relative to the site, [km]
        threshold_pre -> [float,optional,default=[10,5]] Angular distance threshold and slant distance threshold for initial match, [deg,km]
        threshold_deep -> [float,optional,default=[500,1]] Angular distance threshold and slant distance threshold for deep match, [arcsec,km]
        threshold_slope -> [float,optional,default=[1,0.5]] Slope threshold of angular distance and slant distance for deep match, [arcsec/s,km/s]
    Outputs:
        code_match -> [int] match code
        satnum -> [int] ID of matched space object
        disp_match -> [str] description for matched status

    Three types of matching results are summaried as follows    

    | Match Code |   Satnum          | Solution Case      | Status  | What to Do Next    |
    |:----------:|:-----------------:|:------------------:|:-------:|:------------------:|
    |  1         | NORAD ID          | Unique solution    | Success |                    |
    |  0         | None              | No solution        | Failure | increase threshold |
    | -1         | list of NORAD IDs | Multiple solutions | Failure | decrease threshold |    
    """
    # Convert azalt to radec
    az,alt = azalt.T
    a,ecc,inc,raan,argp,nu = orbele_site.T
    x,y,z = spherical_to_cartesian(1,np.deg2rad(alt),np.deg2rad(360-az))
    xyz_RADAR = np.stack([x.value,y.value,z.value]).T
    ECI2RADAR_mat,RADAR2ECI_mat = ECI_RADAR_mat(inc,raan,argp,nu)
    xyz_ECI = Matrix_dot_Vector(RADAR2ECI_mat,xyz_RADAR)
    _r,dec,ra = cartesian_to_spherical(xyz_ECI[:,0],xyz_ECI[:,1],xyz_ECI[:,2])
    radec = np.stack([ra.deg,dec.deg]).T

    # load the objects list in TLE
    sats_list = tle._sats_Satrec

    # Initial matching
    sats_list_filter = match_pre(sats_list,ta,xyz_site,radec,r,threshold_pre)
   
    # Deep matching
    if sats_list_filter:
        sats_list_find = match_deep(sats_list_filter,ta,xyz_site,radec,r,threshold_deep,threshold_slope)
        n = len(sats_list_find)
    else:
        n = 0   
    # Results integration
    if n == 1 :
        code_match = 1 # unique solution
        satnum = sats_list_find[0].satnum
        disp_match = 'Target ID: {:d}'.format(satnum)
    elif n == 0:
        code_match = 0 # no solution
        satnum = None
        disp_match = 'May be a new target! Please try increasing the threshold for further confirmation.'
    else:
        code_match = -1 # mutiple solutions  
        satnum = [sat_find.satnum for sat_find in sats_list_find]
        disp_match = 'Multiple targets associated to the same arc! Please try reducing the threshold for further confirmation.' # The threshold parameters need to be adjusted!

    return code_match,satnum,disp_match   
   
def match_pre(sats_list,ta,xyz_site,radec,r,threshold=[10,10]):
    """
    Initially match the observation arc to space objects in TLE.

    Usage:
        >>> sats_list_filter = match_pre(sats_list,ta,xyz_site,radec,r)
    Inputs:
        sats_list -> objects list in TLE
        ta -> [Astropy Time] time
        xyz_site -> [2D array] cartesian coordinates of site, [km]
        radec -> [2D array] RA and Dec of space object, [deg]
        r -> [array] Slant distance of the space object relative to the site, [km]
        threshold -> [float,optional,default=[10,5]] Angular distance threshold and slant distance threshold for initial match, [deg,km]
    Outputs:
        sats_list_filter -> filtered objects list after matching
    """

    # Extract the data at the middle moment of the arc segment
    ta_mid,xyz_site_mid,radec_mid,r_mid = get_mid_point(ta,xyz_site,radec,r)

    # calculate the cartesian coordinates of space objects in GCRF at the middle moment of the arc
    sats_array = SatrecArray(sats_list)  
    e, xyz_teme, vxyz_teme = sats_array.sgp4(np.array([ta_mid.jd1]),np.array([ta_mid.jd2]))
    gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta_mid)
    xyz_gcrf = (teme2gcrf_mat @ xyz_teme.transpose(1,2,0)).transpose(2,0,1) 

    # Calculate the angular distance between the measured pointing and the calculated pointing of space objects
    site_sat_xyz = (xyz_gcrf - xyz_site_mid).T
    c1 = SkyCoord(site_sat_xyz, unit='km', frame='icrs',representation_type='cartesian')
    c2 = SkyCoord([radec_mid],frame='icrs',unit='deg')
    sep = c1.separation(c2).deg

    # Calculate the difference between the measured slant distance and the calculated slant distance of space objects
    diff_r = np.abs(r_mid - norm(site_sat_xyz,axis=0).flatten())

    # According to the angular distance threshold and slant distance threshold, filter the objects list in TLE
    index, = np.where((sep < threshold[0])&(diff_r < threshold[1]))
    sats_list_filter = [sats_list[i] for i in index]
    
    return sats_list_filter        

def match_deep(sats_list,ta,xyz_site,radec,r,threshold=[2000,5],threshold_slope=[400,2]):
    """
    Deeply match the observation arc to space objects in filtered objects list.

    Usage:
        >>> sats_list_find = match_deep(sats_list,ta,xyz_site,radec,r)
    Inputs:
        sats_list -> filtered objects list after initial matching
        ta -> [Astropy Time] time
        xyz_site -> [2D array] cartesian coordinates of site
        radec -> [2D array] RA and Dec of space object
        r -> [array] Slant distance of the space object relative to the site, [km]
        threshold -> [float,optional,default=[500,1]] Angular distance threshold and slant distance threshold for deep match, [arcsec,km]
        threshold_slope -> [float,optional,default=[1,0.5]] Slope threshold of angular distance and slant distance for deep match, [arcsec/s,km/s]
    Outputs:
        sats_list_find -> objects list after deep matching  
    """

    # calculate the cartesian coordinates of space objects in GCRF at the observation time sequence
    sats_array = SatrecArray(sats_list)  
    e, xyz_teme, vxyz_teme = sats_array.sgp4(ta.jd1,ta.jd2)
    gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta)
    xyz_gcrf = (teme2gcrf_mat @ xyz_teme.transpose(1,2,0)).transpose(2,0,1) 
    site_sat_xyz = (xyz_gcrf - xyz_site).T
    site_sat_x,site_sat_y,site_sat_z = site_sat_xyz

    # Calculate the rms of angular distance
    c1 = SkyCoord(site_sat_x,site_sat_y,site_sat_z, unit='km', frame='icrs',representation_type='cartesian')
    c2 = SkyCoord(radec, frame='icrs',unit='deg')
    sep = c1.separation(c2[:,None]).arcsec
    sep_rms = norm(sep,axis=0)/np.sqrt(len(sep))

    # Calculate the rms of difference between the measured slant distance and the calculated slant distance
    diff_r = r[:,None] - norm(site_sat_xyz,axis=0)
    r_rms = norm(diff_r,axis=0)/np.sqrt(len(r))
    
    sep_flag = sep_rms < threshold[0]
    r_flag = r_rms < threshold[1]

    # Estimate the slope for the angular distance and slant distance
    dta = (ta - ta[0]).sec
    A = np.vstack([dta, np.ones(len(dta))]).T
    slope_sep = lstsq(A, sep,rcond=None)[0][0]
    slope_sep_flag = np.abs(slope_sep) < threshold_slope[0]

    slope_r = lstsq(A, diff_r,rcond=None)[0][0]
    slope_r_flag = np.abs(slope_r) < threshold_slope[1]

    # According to the combined threshold, filter the objects in initial matching
    flag = sep_flag & r_flag & slope_sep_flag & slope_r_flag
    index, = np.where(flag)
    sats_list_find = [sats_list[i] for i in index]

    return sats_list_find