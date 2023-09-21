import numpy as np
from astropy.time import Time
from astropy import units
from astropy.coordinates import spherical_to_cartesian

from ..arcmatch.optical import arcsat_match as arcsat_match_optical
from ..arcmatch.radar import arcsat_match as arcsat_match_radar
from ..utils.preprocessing import lowess_smooth_optical,lowess_smooth_radar
from ..iod.common import getnpoints
from ..transform.frame_trans import ECI_RADAR_mat
from ..utils.math import Matrix_dot_Vector
from .class_iod import IOD

class ArcObs(object):
    """
    class of Observation Arcs(Tracklets)
        Attributes:
            - mode -> [str] Type of the observation arcs. Available options include 'optical' or 'radar', 
            where optical measurements refers to right ascension and declination data, 
            and radar measurements refers to azimuth, altitude, and slant distance data. 
            - t -> [array of str] Time sequence in UTC
            - radec -> [2D array] Ra and Dec of space objects in Site-centered Inertial Reference Frame, [deg]
            - azalt -> [2D array] Az and Alt of space objects in Site-centered RADAR Reference Frame, [deg];
            it is defined as a right-handed reference frame with x-axis pointing Along-track, y-axis Cross-track, and z-axis Radial.
            The azimuth angle is measured clockwise from the x-axis.
            - r -> [array of float] Slant distance of the space object relative to the site, [km]
            - xyz_site -> [3D array] Cartesian coordinates of space objects in GCRF(Geocentric Celestial Reference Frame), [km]
            - _ta_lowess -> [array like] Astropy Time sequence with outliers removed using LOWESS 
            - _radec_lowess -> [2D array] Ra and Dec with outliers removed using LOWESS 
            - _xyz_site_lowess -> [3D array] Cartesian coordinates with outliers removed using LOWESS 
            - _azalt_lowess -> [2D array] Az and Alt with outliers removed using LOWESS 
            - _r_lowess -> [array of float] Slant distance with outliers removed using LOWESS 
            - flag_lowess -> [array of bool] Flag of data points. If False, the data point is an outlier.
        Methods:    
            - lowess_smooth -> Remove outliers in optical/radar data with the method of LOWESS (Locally Weighted Scatterplot Smoothing).
            - arc_match -> Match the observation arc based on optical angle measurement data or radar measurement data(range+angle) to space objects in TLE file.
            - iod -> Extract the necessary information for Initial Orbit Determination(IOD) from optical angle measurement data or radar measurement(range+angle) data.
    """
    def __init__(self,info):
        """
        Create an instance of class ArcObs.

        Usage:
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site})
            >>> # arc_radar = ArcObs({'t':t,'azalt':azalt,'r':r,'xyz_site':xyz_site})
        Inputs:
            info -> [dict] Necessary data for observation arcs in form of dictionary    
        Outputs:
            arc_optical -> An instance of class ArcObs 
        """
        if 'r' in info.keys() and 'azalt' in info.keys():
            info['mode'] = 'radar'
            info['_azalt_lowess'] = info['azalt']
            info['_r_lowess'] = info['r']
            info['_orbele_site_lowess'] = info['orbele_site']
        elif 'radec' in info.keys():
            info['mode'] = 'optical'
            info['_radec_lowess'] = info['radec']
        else:
            raise Exception('Input only supports optical angle measurement data and radar (ranging + angle measurement) data.')    

        info['_ta_lowess'] = Time(info['t'])
        info['_xyz_site_lowess'] = info['xyz_site']

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        return 'Instance of class ArcObs'

    def lowess_smooth(self,frac=0.5):
        """
        Remove outliers in optical/radar data with the method of LOWESS (Locally Weighted Scatterplot Smoothing)
        Here, LOWESS uses a weighted **linear regression** in default.

        Usage:
            >>> arc_optical.lowess_smooth()
        Inputs:
            frac -> [float,optional,default=0.5] The fraction of the data used in local regression. The value of fraction is between 0 and 1.    
        Outputs:
            An instance of class ArcObs with LOWESS processing    
        """
        ta = Time(self.t)
        info = self.info
        if self.mode == 'optical': 
            flag = lowess_smooth_optical(ta,self.radec,frac)
            self._radec_lowess = info['_radec_lowess'] = self.radec[flag]
        elif self.mode == 'radar':
            flag = lowess_smooth_radar(ta,self.azalt,self.r,frac)
            self._r_lowess = info['_r_lowess'] = self.r[flag]
            self._azalt_lowess = info['_azalt_lowess'] = self.azalt[flag]
            self._orbele_site_lowess = info['_orbele_site_lowess'] = self.orbele_site[flag]

        self._ta_lowess = info['_ta_lowess'] = ta[flag]
        self._xyz_site_lowess = info['_xyz_site_lowess'] = self.xyz_site[flag]
        self.flag_lowess = info['flag_lowess'] = flag   

        return None    

    def arc_match(self,tle,threshold_dict=None):
        """
        Match the observation arc based on optical angle measurement data or radar measurement data(range+angle) to space objects in TLE file.

        Usage:
            >>> arc_optical.arc_match(tle)
            >>> arc_radar.arc_match(tle)
        Inputs:
            tle -> An instance of class TLE
            threshold_dict -> [dict,optional,default=None] Threshold for arc matching, such as {'threshold_pre':5,'threshold_deep':200,'threshold_slope':0.5} for optical angle measurement data, 
            {'threshold_pre':[10,5],'threshold_deep':[500,1],'threshold_slope':[1,0.5]} for radar measurement data.
            For case of optical angle measurement data:
                threshold_pre -> [float,optional,default=5] Angular distance threshold for initial match, [deg]
                threshold_deep -> [float,optional,default=200] Angular distance threshold for deep match, [arcsec]
                threshold_slope -> [float,optional,default=0.5] Slope threshold of angular distance for deep match, [arcsec/s]
            For case of radar measurement data:
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
        info = self.info
        ta = self._ta_lowess
        xyz_site = self._xyz_site_lowess

        # For case of optical angle measurement data
        if self.mode == 'optical':
            radec = self._radec_lowess
            if threshold_dict is None:
                threshold_pre = 5
                threshold_deep = 1000 # [arcsec]
                threshold_slope = 200 # [arcsec/s]
            else:
                threshold_pre = threshold_dict['threshold_pre']
                threshold_deep = threshold_dict['threshold_deep']
                threshold_slope = threshold_dict['threshold_slope']
            code_match,satnum,disp_match = arcsat_match_optical(tle,ta,xyz_site,radec,threshold_pre,threshold_deep,threshold_slope)
        
        # For case of radar measurement data(range+angle)
        elif self.mode == 'radar':
            azalt = self._azalt_lowess
            r = self._r_lowess
            orbele_site = self._orbele_site_lowess
            
            if threshold_dict is None:
                threshold_pre = [10,10] # [deg,km]
                threshold_deep = [2000,5] # [arcsec,km]
                threshold_slope = [400,2] # [arcsec/s,km/s]
            else:
                threshold_pre = threshold_dict['threshold_pre']
                threshold_deep = threshold_dict['threshold_deep']
                threshold_slope = threshold_dict['threshold_slope'] 
            code_match,satnum,disp_match = arcsat_match_radar(tle,ta,xyz_site,orbele_site,azalt,r,threshold_pre,threshold_deep,threshold_slope)

        self.code_match = info['code_match'] = code_match
        self.satnum = info['satnum'] = satnum
        self.disp_match = info['disp_match'] = disp_match  

        return None

    def iod(self,body,nd_unit=None):    
        """
        Extract the necessary information for Initial Orbit Determination(IOD) from optical angle-only measurements or radar range+angle measurements.

        Usage:
            >>> # Extract necessary data for IOD from radar range+angle measurements
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test/radar_obs.dat',dtype=str,skiprows=1) # Load the observation file
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> orbele_site = obs_data[:,1:7].astype(float) # Orbital elements(a,ecc,inc,raan,argp,true_anomaly) of the site
            >>> xyz_site = obs_data[:,7:10].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> azalt = obs_data[:,10:12].astype(float) # Azimuth and Altitude angle of space object, [deg]
            >>> r = obs_data[:,12].astype(float) # Slant distance of the space object relative to the site, [km]
            >>> # Load the necessary data
            >>> from orbdtools import ArcObs
            >>> arc_radar = ArcObs({'t':t,'azalt':azalt,'r':r,'xyz_site':xyz_site,'orbele_site':orbele_site}) 
            >>> arc_radar.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_radar.iod(earth)
        Inputs:
            body -> [object of class Body] Central body of attraction. Currently, available options include 'Earth', 'Moon', 'Sun', 'Venus', 'Jupiter', 'Mars' and user-defined celestial bodies.
            nd_unit -> [object of class Body,optional,default=None] Non-dimensional unit system defined by a celestial body. If None, it is the same as that of the central body of attraction.
        Outputs:
            object of class IOD(Initial Orbit Determination)
        """
        # Normalization of units
        if nd_unit is None:
            nd_unit = body
            mu_nd = 1
            Re_nd = 1
        else:  
            mu_nd = body.mu/nd_unit.mu  
            Re_nd = body.Re/nd_unit.Re

        info = self.info
        ta = self._ta_lowess
        xyz_site = self._xyz_site_lowess
        xyz_site_nd = xyz_site / nd_unit._L # Convert to non-dimensional length unit [L_nd]

        # For case of optical angle-only measurements
        if self.mode == 'optical':
            radec = self._radec_lowess
            # Compute the unit vector for the line of sight direction
            losnp = spherical_to_cartesian(1,radec[:,1]*units.deg,radec[:,0]*units.deg)
            losnp = np.stack([los_i.value for los_i in losnp]).T
            # Extract three approximately uniformly distributed data points.
            t3p,xyz_site3p_nd,los3p = getnpoints(3,ta,xyz_site_nd,losnp) 

            # only for the reference vector method
            """
            # Considering that the reference vector method cannot successfully determine the orbit, the code related to this method is temporarily commented out.
            ra_rad = np.deg2rad(radec[:,0])
            alpha_vec = np.array([-np.sin(ra_rad),np.cos(ra_rad),np.zeros(len(ta))]).T
            delta_vec = np.cross(losnp,alpha_vec)
            dict_keys = ['t3p','xyz_site3p_nd','los3p','xyz_site_nd','losnp','_alpha_vec','_delta_vec']
            dict_values = [t3p,xyz_site3p_nd,los3p,xyz_site_nd,losnp,alpha_vec,delta_vec]
            """
            dict_keys = ['t3p','xyz_site3p_nd','los3p','xyz_sitenp_nd','losnp']
            dict_values = [t3p,xyz_site3p_nd,los3p,xyz_site_nd,losnp]

        # For case of radar range+angle measurements
        elif self.mode == 'radar':
            azalt = self._azalt_lowess
            r = self._r_lowess   
            orbele_site = self._orbele_site_lowess

            a,ecc,inc,raan,argp,nu = orbele_site.T # angular variables are in [deg]
            r_nd = r / body._L # Convert slant distance to non-dimensional length unit [L_nd]
            # Compute the unit vector of Line-Of-Sight(LOS)
            losnp = spherical_to_cartesian(1,azalt[:,1]*units.deg,(360-azalt[:,0])*units.deg)
            losnp = np.stack([los_i.value for los_i in losnp]).T
            # Convert the unit vector to ECI from RADAR reference frame
            ECI2RADAR_mat,RADAR2ECI_mat = ECI_RADAR_mat(inc,raan,argp,nu)
            losnp = Matrix_dot_Vector(RADAR2ECI_mat,losnp)
            # Compute the position of the space objects in GCRF(Geocentric Celestial Reference Frame)
            posnp_nd = xyz_site_nd + r_nd[:,None] * losnp

            # Extract three approximately uniformly distributed data points.
            t3p,xyz_site3p_nd,los3p,dis3p_nd,pos3p_nd = getnpoints(3,ta,xyz_site_nd,losnp,r_nd,posnp_nd) 
            dict_keys = ['t3p','xyz_site3p_nd','los3p','dis3p_nd','pos3p_nd','xyz_sitenp_nd','losnp','disnp_nd','posnp_nd']
            dict_values = [t3p,xyz_site3p_nd,los3p,dis3p_nd,pos3p_nd,xyz_site_nd,losnp,r_nd,posnp_nd]

        t1,t2,t3 = t3p
        tof_nd = (t3 - t1).sec/nd_unit._T # Convert time of flight to non-dimensional time unit [T_nd]
        tau_nd = np.array([(t2 - t1).sec,(t3 - t2).sec])/nd_unit._T 

        t0_1 = t3p[1] # Median Epoch at which the orbit elements is estimated. It refers specifically to the median of the observed time series.
        t0_2 = (ta[-1] - ta[0])/2 + ta[0] # Intermediate Epoch at which the orbit elements is estimated. It refers specifically to the middle moment of the observation arc segment.
        t_1_nd = (ta - t0_1).sec / nd_unit._T # Elapsed time since the Median Epoch
        t_2_nd = (ta - t0_2).sec / nd_unit._T # Elapsed time since the Intermediate Epoch

        dict_keys += ['tof_nd','tau_nd','body','nd_unit','mu_nd','Re_nd','t0_1','t0_2','t_1_nd','t_2_nd']
        dict_values += [tof_nd,tau_nd,body,nd_unit,mu_nd,Re_nd,t0_1,t0_2,t_1_nd,t_2_nd]

        add_info = dict(zip(dict_keys, dict_values))
        info.update(add_info)

        return IOD(info)   
