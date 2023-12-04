import numpy as np

class COD(object):
    """
    Class of COD(Cataloging Orbit Determination).
        Attributes:
            - mode -> [str] Type of the observations. Available options include 'optical' and 'radar'. 
            Optical measurements refer specifically to right ascension(RA) and declination(DEC) data, 
            and radar measurements refer specifically to azimuth, altitude, and slant distance data. 
            - t -> [array of str] Time sequence in UTC
            - tof -> [float] Time of flight in [sec]
            - xyz_site -> [3D array] Cartesian coordinates of the site in GCRF(Geocentric Celestial Reference Frame), [km]
            - losnp -> [2D array with shape of nx3] Line-Of-Sight(LOS) vector of the space object relative to the site
            - flag_lowess -> [array of bool] Validity flag of data points. If False, the data point is regarded as an outlier.
            - radec -> [2D array] Ra and Dec of space objects in Site-centered Inertial Reference Frame, with unit of [deg]
            - azalt -> [2D array] Az and Alt of space objects in Site-centered RADAR Reference Frame, with unit of [deg].
            The RADAR reference frame is defined as a right-handed reference frame with x-axis pointing Along-track, y-axis Cross-track, and z-axis Radial.
            The azimuth angle is measured clockwise from the x-axis.
            - r -> [array of float] Slant distance of the space object relative to the site, with unit of [km]
            - posnp -> [2D array] Cartesian coordinates of space objects in GCRF(Geocentric Celestial Reference Frame), [km]
            - orbele_site -> [2D array] Orbital elements(a, e, i, Ω, ω, ν) of the site
            - _ta_lowess -> [array like, object of Astropy Time] Time sequence with outliers removed using LOWESS 
            - _xyz_site_lowess -> [2D array of float] Cartesian coordinates of site with outliers removed using LOWESS 
            - _losnp_lowess -> [2D array] Line of Sight(LOS) in GCRF(Geocentric Celestial Reference Frame) with outliers removed using LOWESS 
            - _radec_lowess -> [2D array of float] Ra and Dec with outliers removed using LOWESS 
            - _azalt_lowess -> [2D array] Az and Alt with outliers removed using LOWESS 
            - _r_lowess -> [array of float] Slant distance with outliers removed using LOWESS 
            - _posnp_lowess -> [2D array] Cartesian coordinates of space objects in GCRF(Geocentric Celestial Reference Frame) with outliers removed using LOWESS 
            - _orbele_site_lowess -> [2D array] Orbital elements(a, e, i, Ω, ω, ν) of the site with outliers removed using LOWESS
            - df -> [Dataframe of Pandas] Dataframe of improved orbital elements with keys and values as follows
                epoch -> [str] Epoch of the orbital elements in UTC
                a -> [float] Semimajor in [Length unit], where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km. 
                Considering that the orbit elements are only initial values, [Length unit] can take slightly different values, such as the equatorial radius defined by the WGS84 model.
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [deg]
                raan -> [float] Longitude of ascending node, [deg]
                argp -> [float] Argument of perigee, [deg]
                nu -> [float] True anomaly, [deg]. 
                M -> [float] Mean anomaly, [deg]. 
                h -> [float] Magnitude of angular momentum
                bstar -> [float ] The drag term, or radiation pressure coefficient in unit of [1/earth radii]. Its absolute value is usually less than 1.
                status -> [str] Convergence state of the estimation: 'success' or 'failed'
            - method -> [str] Method for COD. Available options include 'SGP4' and 'DSST'. 
            - rms -> [float] RMS of O-C   
            - ref -> [str] Reference frame in which the improved orbital elements is defined. For SGP4, it should be TEME by default; otherwise, it should be one of GCRF/J2000/ECI by default.
            - satid -> [int] Cataloging ID of the space object  
            - tle_str -> [list of str] Two line elements 
    """
    def __init__(self,info):
        """
        Initialize an instance of class COD.
        """
        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a more information-rich string representation of the COD object.
        """
        return '<COD object: {:s} Start Time = {:s}Z TOF ≈ {:.0f}s Method = {:s} ref = {:s} Elements = {}>'.format(self.mode.upper(),self._ta_lowess[0].isot,self.tof,self.method,self.ref,self.df.to_dict('records'))