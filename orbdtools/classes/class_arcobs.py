import numpy as np
from astropy.time import Time
from astropy import units
from astropy.coordinates import spherical_to_cartesian
from copy import deepcopy
import pandas as pd

from ..arcmatch.optical import arcsat_match as arcsat_match_optical
from ..arcmatch.radar import arcsat_match as arcsat_match_radar
from ..utils.preprocessing import lowess_smooth_optical,lowess_smooth_radar
from ..iod.common import getnpoints
from ..transform.frame_trans import ECI_RADAR_mat
from ..utils.math import Matrix_dot_Vector
from ..utils import Const
from ..cod.sgp4.sgp4_od.optical import sgp4_od_optical
from ..cod.sgp4.sgp4_od.radar import sgp4_od_radar
from .class_iod import IOD
from .class_cod import COD

class ArcObsArray(object):
    """
    Class of Observation Arc(Tracklet) Array
        Attributes:
            - arcs -> Array of ArcObs object
            - mode -> Observation type of the arcs
        Methods:
            - cod_sgp4 -> Cataloging orbit determination using combined optical angle measurements and radar range+angle measurements with SGP4 propagator. 
    """
    def __init__(self,info):
        """
        Initialize an instance of class ArcObsArray.
        """
        for key in info.keys():
            setattr(self, key, np.array(info[key]))

    def __repr__(self):
        """
        Returns a more information-rich string representation of the ArcObsArray object.
        """
        n = len(self.arcs)
        arcs_mode = [mode.upper() for mode in self.mode]
        return '<ArcObsArray object: LENGTH = {:d} OBSTYPE = {}>'.format(n,arcs_mode)        

class ArcObs(object):
    """
    class of Observation Arcs(Tracklets)
        Attributes:
            - mode -> [str] Type of the observation arcs. Available options include 'optical' or 'radar', 
            where optical measurements refers to right ascension and declination data, 
            and radar measurements refers to azimuth, altitude, and slant distance data. 
            - t -> [array of str] Time sequence in UTC
            - tof -> [float] Time of flight in [sec]
            - xyz_site -> [2D array] Cartesian coordinates of the site in GCRF(Geocentric Celestial Reference Frame), [km]
            - losnp -> [2D array] Line of Sight(LOS) in GCRF(Geocentric Celestial Reference Frame)
            - flag_lowess -> [array of bool] Validity flag of data points. If False, the data point is an outlier.
            - radec -> [2D array] Ra and Dec of space objects in Site-centered Inertial Reference Frame, [deg]
            - azalt -> [2D array] Az and Alt of space objects in Site-centered RADAR Reference Frame, [deg];
            it is defined as a right-handed reference frame with x-axis pointing Along-track, y-axis Cross-track, and z-axis Radial.
            The azimuth angle is measured clockwise from the x-axis.
            - r -> [array of float] Slant distance of the space object relative to the site, [km]
            - posnp -> [2D array] Cartesian coordinates of space objects in GCRF(Geocentric Celestial Reference Frame), [km]
            - orbele_site -> [2D array] Orbital elements(a, e, i, Ω, ω, ν) of the site
            - _ta_lowess -> [array like] Astropy Time sequence with outliers removed using LOWESS 
            - _xyz_site_lowess -> [3D array] Cartesian coordinates of the site with outliers removed using LOWESS 
            - _losnp_lowess -> [2D array] Line of Sight(LOS) in GCRF(Geocentric Celestial Reference Frame) with outliers removed using LOWESS 
            - _radec_lowess -> [2D array] Ra and Dec with outliers removed using LOWESS 
            - _azalt_lowess -> [2D array] Az and Alt with outliers removed using LOWESS 
            - _r_lowess -> [array of float] Slant distance with outliers removed using LOWESS
            - _posnp_lowess -> [2D array] Cartesian coordinates of space objects in GCRF(Geocentric Celestial Reference Frame) with outliers removed using LOWESS 
            - _orbele_site_lowess -> [2D array] Orbital elements(a, e, i, Ω, ω, ν) of the site with outliers removed using LOWESS
        Methods:    
            - lowess_smooth -> Remove outliers in optical/radar data with the method of LOWESS (Locally Weighted Scatterplot Smoothing).
            - arc_match -> Match the observation arc based on optical angle measurement data or radar measurement data(range+angle) to space objects in TLE file.
            - iod -> Extract the necessary information for Initial Orbit Determination(IOD) from optical angle measurement data or radar measurement(range+angle) data.
            - cod_sgp4 -> Cataloging orbit determination using optical angle measurements or radar range+angle measurements with SGP4 propagator. 
    """
    def __init__(self,info):
        """
        Create an instance of class ArcObs.

        Usage:
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site})
            >>> # arc_radar = ArcObs({'t':t,'azalt':azalt,'r':r,'xyz_site':xyz_site,'orbele_site':orbele_site})
        Inputs:
            info -> [dict] Necessary information for initializing observation arcs in form of dictionary    
        Outputs:
            arc -> An instance of class ArcObs 
        """
        info['_ta_lowess'] = _ta_lowess = Time(info['t'])
        xyz_site = info['_xyz_site_lowess'] = info['xyz_site']
        info['flag_lowess'] = np.ones_like(info['t']).astype(bool)
        tof = (_ta_lowess[-1] - _ta_lowess[0]).sec
        info['tof'] = tof

        if 'r' in info.keys() and 'azalt' in info.keys():
            info['mode'] = 'radar'
            azalt = info['_azalt_lowess'] = info['azalt']
            r = info['_r_lowess'] = info['r']
            orbele_site = info['_orbele_site_lowess'] = info['orbele_site']

            # Compute the unit vector of Line-Of-Sight(LOS)
            losnp = spherical_to_cartesian(1,azalt[:,1]*units.deg,(360-azalt[:,0])*units.deg)
            losnp = np.stack([los_i.value for los_i in losnp]).T
            a,ecc,inc,raan,argp,nu = orbele_site.T # angular variables are in [deg]

            # Convert the unit vector to ECI from RADAR reference frame
            ECI2RADAR_mat,RADAR2ECI_mat = ECI_RADAR_mat(inc,raan,argp,nu)
            losnp = Matrix_dot_Vector(RADAR2ECI_mat,losnp)

            # Compute the position of the space objects in GCRF(Geocentric Celestial Reference Frame)
            info['_posnp_lowess'] = info['posnp'] = posnp = xyz_site + r[:,None] * losnp
        elif 'radec' in info.keys():
            info['mode'] = 'optical'
            radec = info['_radec_lowess'] = info['radec']

            # Compute the unit vector for the line of sight direction
            losnp = spherical_to_cartesian(1,radec[:,1]*units.deg,radec[:,0]*units.deg)
            losnp = np.stack([los_i.value for los_i in losnp]).T
        else:
            raise Exception('Only optical angle measurement data or radar(ranging + angle) measurement data are supported.')

        info['_losnp_lowess'] = info['losnp'] = losnp

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a more information-rich string representation of the ArcObs object.
        """
        n_data = len(self.t)
        n_outlier = np.invert(self.flag_lowess).sum()
        return '<ArcObs object: {:s} Start Time = {:s}Z TOF ≈ {:.0f}s N_DATA = {:d} N_OUTLIER = {:d}>'.format(self.mode.upper(),self._ta_lowess[0].isot,self.tof,n_data,n_outlier)

    def lowess_smooth(self,frac=0.5):
        """
        Remove outliers in optical/radar data with the method of LOWESS (Locally Weighted Scatterplot Smoothing)
        Here, LOWESS uses a weighted **linear regression** in default.

        Usage:
            >>> arc_optical.lowess_smooth()
            >>> # arc_radar.lowess_smooth()
        Inputs:
            frac -> [float,optional,default=0.5] The fraction of the data used in local regression. The value should be between 0 and 1.    
        Outputs:
            An instance of class ArcObs with LOWESS processing    
        """
        ta = Time(self.t)
        if self.mode == 'optical': 
            flag = lowess_smooth_optical(ta,self.radec,frac)
            self._radec_lowess = self.radec[flag]
        elif self.mode == 'radar':
            flag = lowess_smooth_radar(ta,self.azalt,self.r,frac)
            self._r_lowess = self.r[flag]
            self._azalt_lowess = self.azalt[flag]
            self._orbele_site_lowess = self.orbele_site[flag]
            self._posnp_lowess = self.posnp[flag]

        self._ta_lowess = _ta_lowess = ta[flag]
        self._xyz_site_lowess = self.xyz_site[flag]
        self._losnp_lowess = self.losnp[flag]
        self.flag_lowess = flag   
        self.tof = (_ta_lowess[-1] - _ta_lowess[0]).sec

        return None  

    def fuse(self,arcs):
        """
        Fuse multiple tracklets of different type, and make a sort according to the starting time of the arcs.

        Usage:
            >>> arcs_array = arc_optical1.fuse([arc_optical2,arc_radar1])
        Inputs:
            arcs -> [list of ArcObs objects] List of multiple tracklets    
        Outputs:
            arcs_array -> An instance of class ArcObsArray   
        """
        arcs = np.append(self,arcs)
        arcs_t = [arc.t[0] for arc in arcs]
        arcs_index = np.argsort(arcs_t)
        arcs = arcs[arcs_index]  
        mode = [arc.mode for arc in arcs]
        info = {'arcs':arcs,'mode':mode}
        return ArcObsArray(info) 

    def join(self,arcs):
        """
        Simply join multiple tracklets of same type into one long arc.

        Usage:
            >>> arc_optical = arc_optical1.join([arc_optical2,arc_optical2])
        Inputs:
            arcs -> [list of ArcObs objects] List of multiple tracklets    
        Outputs:
            arc_optical -> An instance of class ArcObs   
        """
        arcs = self.fuse(arcs)
        arcs_list = arcs.arcs
        arc0 = deepcopy(arcs_list[0])

        if len(set(arcs.mode)) != 1:
            raise Exception("Cannot join multiple tracklets of different type into one long arc")

        arcs_mode = arc0.mode
        
        if arcs_mode == 'optical':
            for arc in arcs_list[1:]:
                arc0.t = np.append(arc0.t,arc.t)
                arc0.xyz_site = np.append(arc0.xyz_site,arc.xyz_site,axis=0)
                arc0._xyz_site_lowess = np.append(arc0._xyz_site_lowess,arc._xyz_site_lowess,axis=0)
                arc0._ta_lowess = np.append(arc0._ta_lowess,arc._ta_lowess)
                arc0.flag_lowess = np.append(arc0.flag_lowess,arc.flag_lowess)
                arc0.radec = np.append(arc0.radec,arc.radec,axis=0)
                arc0._radec_lowess = np.append(arc0._radec_lowess,arc._radec_lowess,axis=0)  
                arc0.losnp = np.append(arc0.losnp,arc.losnp,axis=0)
                arc0._losnp_lowess = np.append(arc0._losnp_lowess,arc._losnp_lowess,axis=0)
        elif arcs_mode == 'radar':  
            for arc in arcs_list[1:]:
                arc0.t = np.append(arc0.t,arc.t)
                arc0.xyz_site = np.append(arc0.xyz_site,arc.xyz_site,axis=0)
                arc0._xyz_site_lowess = np.append(arc0._xyz_site_lowess,arc._xyz_site_lowess,axis=0)
                arc0._ta_lowess = np.append(arc0._ta_lowess,arc._ta_lowess)
                arc0.flag_lowess = np.append(arc0.flag_lowess,arc.flag_lowess)  
                arc0.azalt = np.append(arc0.azalt,arc.azalt,axis=0)
                arc0._azalt_lowess = np.append(arc0._azalt_lowess,arc._azalt_lowess,axis=0)
                arc0.r = np.append(arc0.r,arc.r)
                arc0._r_lowess = np.append(arc0._r_lowess,arc._r_lowess)
                arc0.losnp = np.append(arc0.losnp,arc.losnp,axis=0)
                arc0._losnp_lowess = np.append(arc0._losnp_lowess,arc._losnp_lowess,axis=0) 
                arc0.posnp = np.append(arc0.posnp,arc.posnp,axis=0)
                arc0._posnp_lowess = np.append(arc0._posnp_lowess,arc._posnp_lowess,axis=0) 
                arc0.orbele_site = np.append(arc0.orbele_site,arc.orbele_site,axis=0)
                arc0._orbele_site_lowess = np.append(arc0._orbele_site_lowess,arc._orbele_site_lowess,axis=0)

        arc0._ta_lowess = Time(arc0._ta_lowess)    
        arc0.tof = (arc0._ta_lowess[-1] - arc0._ta_lowess[0]).sec
        
        return arc0        

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
            code_match -> [int] Match code
            satid -> [int] ID of matched space object
            disp_match -> [str] Description for matched status

        Three types of matching results are summaried as follows    

        | Match Code |   Satid           | Solution Case      | Status  | What to Do Next    |
        |:----------:|:-----------------:|:------------------:|:-------:|:------------------:|
        |  1         | NORAD ID          | Unique solution    | Success |                    |
        |  0         | None              | No solution        | Failure | increase threshold |
        | -1         | list of NORAD IDs | Multiple solutions | Failure | decrease threshold |    
        """
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
            code_match,satid,disp_match = arcsat_match_optical(tle,ta,xyz_site,radec,threshold_pre,threshold_deep,threshold_slope)
        
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
            code_match,satid,disp_match = arcsat_match_radar(tle,ta,xyz_site,orbele_site,azalt,r,threshold_pre,threshold_deep,threshold_slope)

        self.code_match = code_match
        self.satid = satid
        self.disp_match = disp_match  

        return None

    def iod(self,body,nd_unit=None):    
        """
        Extract the necessary information for Initial Orbit Determination(IOD) from optical angle-only measurements or radar range+angle measurements.

        Usage:
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test/radar_obs.dat',dtype=str,skiprows=1) # Load the observation file
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> orbele_site = obs_data[:,1:7].astype(float) # Orbital elements(a,ecc,inc,raan,argp,true_anomaly) of the site
            >>> xyz_site = obs_data[:,7:10].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> azalt = obs_data[:,10:12].astype(float) # Azimuth and Altitude angle of space object, [deg]
            >>> r = obs_data[:,12].astype(float) # Slant distance of the space object relative to the site, [km]
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
            arc_iod -> Object of class IOD(Initial Orbit Determination)
        """
        # Normalization of units
        if nd_unit is None:
            nd_unit = body
            mu_nd = 1
            Re_nd = 1
        else:  
            mu_nd = body.mu/nd_unit.mu  
            Re_nd = body.Re/nd_unit.Re

        info = self.__dict__
        ta = self._ta_lowess
        xyz_site = self._xyz_site_lowess
        xyz_site_nd = xyz_site / nd_unit._L # Convert to non-dimensional length unit [L_nd]
        losnp = self.losnp

        # For case of optical angle-only measurements
        if self.mode == 'optical':
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
            dict_keys = ['t3p','xyz_site3p_nd','los3p','xyz_sitenp_nd']
            dict_values = [t3p,xyz_site3p_nd,los3p,xyz_site_nd]

        # For case of radar range+angle measurements
        elif self.mode == 'radar':
            r_nd = self._r_lowess / body._L # Convert slant distance to non-dimensional length unit [L_nd]
            posnp_nd = self._posnp_lowess / body._L

            # Extract three approximately uniformly distributed data points.
            t3p,xyz_site3p_nd,los3p,r3p_nd,pos3p_nd = getnpoints(3,ta,xyz_site_nd,losnp,r_nd,posnp_nd) 
            dict_keys = ['t3p','xyz_site3p_nd','los3p','r3p_nd','pos3p_nd','xyz_sitenp_nd','losnp','r_nd','posnp_nd']
            dict_values = [t3p,xyz_site3p_nd,los3p,r3p_nd,pos3p_nd,xyz_site_nd,losnp,r_nd,posnp_nd]

        t1,t2,t3 = t3p
        tof_nd = self.tof/nd_unit._T # Convert time of flight to non-dimensional time unit [T_nd]
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

    def cod_sgp4(self,ele0_dict=None,tle=None,satid=0,**kwargs):    
        """
        Cataloging orbit determination using optical angle measurements or radar range+angle measurements with SGP4 propagator. 

        Usage:    
            >>> satid = 22694
            >>> # For the case that initial orbital elements is known
            >>> ele0_dict = {'epoch': '2022-03-24T19:43:11.000Z','a': 6.85997,'ecc': 0.02590,'inc': 18.53046,'raan': 3.51368,'argp': 170.27359,'M': 14.93570}
            >>> arc_cod = arc_optical.cod_sgp4(ele0_dict,satid)
            >>> # For the case that a TLE database is provided
            >>> from orbdtools import TLE
            >>> tle_file = 'test6/tle_20220722.txt'
            >>> tle = TLE.from_file(tle_file)
            >>> arc_cod = arc_optical.cod_sgp4(tle=tle)
            >>> print(arc_cod.df)
            >>> print(arc_cod.tle_str)
            >>> print(arc_cod.rms)
        Inputs:
            ele0_dict -> [dict,optional,default=None] Initial orbit elements in form of dictionary, where at least the following keys need to be included: 
                epoch -> [str or Astropy Time object] Epoch of the orbital elements in UTC
                a -> [float] Semimajor in [Length unit], where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km. 
                Considering that the orbit elements are only initial values, [Length unit] can take slightly different values, such as the equatorial radius defined by the WGS84 model.
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [deg]
                raan -> [float] Longitude of ascending node, [deg]
                argp -> [float] Argument of perigee, [deg]
                M -> [float] Mean anomaly, [deg]
                If ele0_dict is None, the TLE database should be provided for retrieving initial orbit elements.
            tle -> [TLE object,optional,default=None] The TLE database. If None, ele0_dict should be provided.
            satid -> [int,optional,default=0] Satellite Catalog number, such as NORAD ID.
            kwargs -> Extended parameters for generating TLE:
                reff -> [str,optional,default='GCRF'] Reference Frame at which the initial orbit elements is defined. Available options include 'GCRF','J2000', 'ECI', and 'TEME'.
                Here, 'GCRF', 'J2000', and 'ECI' are treated as equivalent.
                bstar -> [float,optional,default=0] The drag term, or radiation pressure coefficient in unit of [1/earth radii]. Its absolute value is usually less than 1.
                nddot -> [float,optional,default=0] One sixth the second derivative of mean motion in unit of [rad/min^3]. Note that this item is only used for SGP, not SGP4. 
                classification -> [str,optional,default='U'] Classification (U: unclassified, C: classified, S: secret)
                intldesg -> [str,optional,default='YYXXXA'] International Designator, where YY is the last two digits of launch year, XXX is the launch number of the year, and A is the piece of the launch.
                elnum -> [int,optional,default=0] Element set number. Incremented when a new TLE is generated for this object.
                revnum -> [int,optional,default=1000] Revolution number at epoch 
        Outputs:
            arc_cod -> Object of class COD      
        """
                
        if tle is not None:
            if satid not in tle.df['noradid'].to_list(): 
                raise Exception('The space object with NORADID of {:d} is not in the TLE database!'.format(satid))  
            ta0,ele0,params = self._tle_update(tle,satid)
        else:    
            if ele0_dict is not None:
                ta0,ele0,params = self._tle_generate(ele0_dict,satid,**kwargs)
            else:    
                raise Exception("Neither Initial orbital elements nor TLE data is provided.")

        ta = self._ta_lowess
        xyz_site = self._xyz_site_lowess
        losnp = self._losnp_lowess
        
        if self.mode == 'optical':
            ele_dict,tle_str,rms = sgp4_od_optical(ele0,ta0,ta,xyz_site,losnp,params)
        elif self.mode == 'radar':
            posnp = self._posnp_lowess
            ele_dict,tle_str,rms = sgp4_od_radar(ele0,ta0,ta,xyz_site,posnp,params)

        ele_df = pd.DataFrame([ele_dict])

        info = self.__dict__.copy()
        info['method'] = 'SGP4'
        info['ref'] = 'TEME'
        info['satid'] = satid
        info['tle_str'] = tle_str
        info['df'] = ele_df
        info['rms'] = rms
                
        return COD(info)
  
    def _tle_update(self,tle,satid):    
        """
        Update the TLE using optical angle measurements or radar range+angle measurements with SGP4 propagator. 

        Usage:
            >>> from orbdtools import TLE
            >>> tle_file = 'test6/tle_20220722.txt'
            >>> tle = TLE.from_file(tle_file)
            >>> satid = 22694
            >>> ta0,ele0,params = arc_optical._tle_update(tle,satid)
        Inputs:
            tle -> [TLE object,optional,default=None] The TLE database. If None, ele0_dict should be provided.
            satid -> [int,optional,default=0] Satellite Catalog number, such as NORAD ID.
        Outputs:
            ta0 -> [object of Astropy Time] The last epoch moment of the observation  
            ele0 -> [array of float]  Mean orbit elements at ta0 by propagating the given TLE
            params -> [llist] Extended parameters for updating TLE
        """    
        if satid not in tle.df['noradid'].to_list(): 
            raise Exception('The space object with NORAD ID of {:d} is not in the TLE database!'.format(satid))    

        tle_r = tle.retrieve([satid])
        ta1 = Time(tle_r.df['epoch'].values.astype(str))[0]

        ta0 = self._ta_lowess[-1]
        tle_epoch = tle_r.atEpoch(ta0)
        ele0 = tle_epoch.df.loc[:,'a':'M'].values[0] 
        
        reff = 'TEME'
        sat_model = tle_epoch._sats_EarthSatellite[0].model    
        bstar = sat_model.bstar   
        nddot = sat_model.nddot
        classification = sat_model.classification
        intldesg = sat_model.intldesg
        elnum = sat_model.elnum + 1
        revnum = sat_model.revnum + int((ta0 - ta1).sec/60 * sat_model.no/Const.twopi)
        params = [satid,reff,bstar,nddot,classification,intldesg,elnum,revnum]

        return ta0,ele0,params

    def _tle_generate(self,ele0_dict,satid,**kwargs):    
        """
        Generating the TLE using optical angle measurements or radar range+angle measurements with SGP4 propagator. 

        Usage:    
            >>> ele0_dict = {'epoch': '2022-03-24T19:43:11.000Z','a': 6.85997,'ecc': 0.02590,'inc': 18.53046,'raan': 3.51368,'argp': 170.27359,'M': 14.93570}
            >>> satid = 22694
            >>> ta0,ele0,params = arc_optical.cod_sgp4(ele0_dict,satid)
        Inputs:
            ele0_dict -> [dict,optional,default=None] Initial orbit elements in form of dictionary, where at least the following keys need to be included: 
                epoch -> [str or Astropy Time object] Epoch of the orbital elements in UTC
                a -> [float] Semimajor in [Length unit], where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km. 
                Considering that the orbit elements are only initial values, [Length unit] can take slightly different values, such as the equatorial radius defined by the WGS84 model.
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [deg]
                raan -> [float] Longitude of ascending node, [deg]
                argp -> [float] Argument of perigee, [deg]
                M -> [float] Mean anomaly, [deg]
                If ele0_dict is None, the TLE database should be provided for retrieving initial orbit elements.
            satid -> [int,optional,default=0] Satellite Catalog number, such as NORAD ID.
            kwargs -> Extended parameters for generating TLE:
                reff -> [str,optional,default='GCRF'] Reference Frame at which the initial orbit elements is defined. Available options include 'GCRF','J2000', 'ECI', and 'TEME'.
                Here, 'GCRF', 'J2000', and 'ECI' are treated as equivalent.
                bstar -> [float,optional,default=0] The drag term, or radiation pressure coefficient in unit of [1/earth radii]. Its absolute value is usually less than 1.
                ndot -> [float,optional,default=None] Ballistic coefficient or one half the first time derivative of the mean motion in unit of [rad/min^2]. If None, it will be calculated by an approximate formula. Note that this item is only used for SGP, not SGP4. 
                nddot -> [float,optional,default=0] One sixth the second derivative of mean motion in unit of [rad/min^3]. Note that this item is only used for SGP, not SGP4. 
                classification -> [str,optional,default='U'] Classification (U: unclassified, C: classified, S: secret)
                intldesg -> [str,optional,default='YYXXXA'] International Designator, where YY is the last two digits of launch year, XXX is the launch number of the year, and A is the piece of the launch.
                elnum -> [int,optional,default=0] Element set number. Incremented when a new TLE is generated for this object.
                revnum -> [int,optional,default=1000] Revolution number at epoch 
            Outputs:
                ta0 -> [object of Astropy Time] Epoch of the initial orbit element
                ele0 -> [array of float]  Initial orbit element
                params -> [llist] Extended parameters for generating TLE    
        """
        ta0 = Time(ele0_dict['epoch'])
        ele0 = [ele0_dict[key] for key in ['a','ecc','inc','raan','argp','M']]

        if 'reff' in kwargs: 
            reff = kwargs['reff']  
        else:
            reff = 'GCRF'

        if 'bstar' in kwargs:  
            bstar = kwargs['bstar']
        else:    
            bstar = 0   

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
            revnum = 1000

        params = [satid,reff,bstar,nddot,classification,intldesg,elnum,revnum]

        return ta0,ele0,params           
