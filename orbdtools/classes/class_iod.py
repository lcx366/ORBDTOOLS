from ..iod.radar.gibbs import gibbs_estimate
from ..iod.radar.multifitting import ellipse_fitting
from ..iod.angular.circular import solve_circular
from ..iod.angular.gauss import gauss_estimate
from ..iod.angular.laplace import laplace_estimate
from ..iod.angular.doubleR import doubleR_estimate
from ..iod.angular.gooding import gooding_estimate
from ..iod.angular.multilaplace import multilaplace_estimate
from ..iod.angular.fg_series import fg_series_optical
from ..iod.radar.fg_series import fg_series_radar
from ..iod.common import to_ele_dict_radar,to_ele_dict_optical

from ..transform.kep_rv_trans import coe2rv

#from ..iod.angular.ref_vec import ref_vec_estimate
#from ..iod.angular.karimi import karimi_estimate

class IOD(object):
    """
    Class of IOD(Initial Orbit Determination).
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
            - _xyz_site_lowess -> [2D array of float] Cartesian coordinates of the site with outliers removed using LOWESS 
            - _losnp_lowess -> [2D array] Line of Sight(LOS) in GCRF(Geocentric Celestial Reference Frame) with outliers removed using LOWESS 
            - _radec_lowess -> [2D array of float] Ra and Dec with outliers removed using LOWESS 
            - _azalt_lowess -> [2D array] Az and Alt with outliers removed using LOWESS 
            - _r_lowess -> [array of float] Slant distance with outliers removed using LOWESS 
            - _posnp_lowess -> [2D array] Cartesian coordinates of space objects in GCRF(Geocentric Celestial Reference Frame) with outliers removed using LOWESS 
            - _orbele_site_lowess -> [2D array] Orbital elements(a, e, i, Ω, ω, ν) of the site with outliers removed using LOWESS
            - tof_nd -> [float] Time of flight in unit of [T_nd]
            - xyz_sitenp_nd -> [2D array with shape of nx3] Cartesian coordinates of the site in unit of [L_nd]
            - r_nd -> [array-like] Slant distance of the space object relative to the site in unit of [L_nd]
            - posnp_nd -> [2D array with shape of nx3] Cartesian coordinates of the space object in GCRF in unit of [L_nd]
            - t3p -> [array with 3 Astropy Time object] The start, midpoint, and end times of the observation
            - xyz_site3p_nd -> [2D array with shape of 3x3] Cartesian coordinates of the site at t3p in GCRF in unit of [L_nd]
            - los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at t3p
            - r3p_nd -> [array with 3 elements] Slant distance of the space object relative to the site at t3p in unit of [L_nd]
            - pos3p_nd -> [2D array with 3 shape of 3x3] Cartesian coordinates of the space object in GCRF at t3p in unit of [L_nd]
            - tau_nd -> [tuple of float] (t2-t1,t3-t2) in unit of [T_nd], where (t1,t2,t3) represents the start, midpoint, and end times of the observation
            - body -> [object of class Body] Central body of attraction
            - nd_unit -> [object of class Body] Non-dimensional unit system defined by a celestial body
            - mu_nd -> [float] GM of the central body in non-dimensional unit
            - Re_nd -> [float] Equatorial radius of the central body in unit of [L_nd]
            - t0_1 -> [Astropy Time object] Median Epoch at which the orbit elements is estimated. It refers specifically to the median of the observed time series.
            - t0_2 -> [Astropy Time object] Intermediate Epoch at which the orbit elements is estimated. It refers specifically to the middle moment of the observation arc segment.
            - t_1_nd -> [array of float] Elapsed time since the Median Epoch
            - t_2_nd -> [array of float] Elapsed time since the Intermediate Epoch
            - df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of ascending node, [rad] or [deg]
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
            - rms -> [float] RMS of O-C    
            - method -> [str] Method for IOD    
        Methods:    
            - gibbs -> Estimate the classical orbital elements from three-points radar range+angle measurements using Gibbs/Herrick-Gibbs Method.
            - ellipse -> Estimate the classical orbital elements from multiple-points radar range+angle measurements using Elliptical Orbit Fitting Method.
            - circular -> Estimate the classical orbital elements from optical angle-only measurements using Near-Circular Orbit Hypothesis method.
            - gauss -> Estimate the classical orbital elements from optical angle-only measurements using Gauss method.
            - laplace -> Estimate the classical orbital elements from optical angle-only measurements using Laplace method.
            - multilaplace -> Estimate the classical orbital elements from optical angle-only measurements using Multiple-points Laplace method.
            - doubleR -> Estimate the classical orbital elements from optical angle-only measurements using Double-R method.
            - gooding -> Estimate the classical orbital elements from optical angle-only measurements using Gooding method.
            - fg_series -> Estimate the classical orbital elements from optical angle-only measurements or radar range+angle measurements using FG-Series method.
        Note: Attributes with '_nd' suffix mean non-dimensional or dimensionless. The non-dimensional unit of GM, length, and time are represented by [mu_nd], [L_nd], [T_nd] respectively. For more information, please refer to orbdtools/utils/Const.py    
    """
    def __init__(self,info):
        """
        Initialize an instance of class IOD.
        """
        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a more information-rich string representation of the IOD object.
        """
        if hasattr(self, "method"):
            return '<IOD object: {:s} Start Time = {:s}Z TOF ≈ {:.0f}s Central Attraction = {} Method = {:s} Elements = {}>'.format(self.mode.upper(),self._ta_lowess[0].isot,self.tof,self.body,self.method,self.df.to_dict('records'))
        else:
            return '<IOD object: {:s} Start Time = {:s}Z TOF ≈ {:.0f}s Central Attraction = {}>'.format(self.mode.upper(),self._ta_lowess[0].isot,self.tof,self.body)    

    def gibbs(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from radar range+angle measurements using the Gibbs/Herrick-Gibbs method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test/radar_obs.dat',dtype=str,skiprows=1) 
            >>> # Extract necessary data for IOD from radar range+angle measurements
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> orbele_site = obs_data[:,1:7].astype(float) # Orbital elements(a,ecc,inc,raan,argp,true_anomaly) of the site
            >>> xyz_site = obs_data[:,7:10].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> azalt = obs_data[:,10:12].astype(float) # Azimuth and Altitude angle of space object, [deg]
            >>> r = obs_data[:,12].astype(float) # Slant distance of the space object relative to the site, [km]
            >>> # Load the necessary data to ArcObs
            >>> from orbdtools import ArcObs
            >>> arc_radar = ArcObs({'t':t,'azalt':azalt,'r':r,'xyz_site':xyz_site,'orbele_site':orbele_site}) 
            >>> arc_radar.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_radar.iod(earth)
            >>> arc_iod.gibbs(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of ascending node, [rad] or [deg]
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references:
            Vallado D. Fundamentals of Astrodynamics and Applications(4th)[M], Microcosm Press, 2013.  
            Curtis H D. Orbital Mechanics for Engineering Students: Revised 4th edition[M]. Butterworth-Heinemann, 2020.       
        """
        if self.mode == 'optical':
            raise Exception('The method of Gibbs/Herrick-Gibbs is not applicable to optical angle-only measurements for initial orbit determination.')

        pos3p_nd = self.pos3p_nd
        posnp_nd = self.posnp_nd
        tof_nd = self.tof_nd
        tau_nd = self.tau_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd

        ele = gibbs_estimate(mu_nd,tof_nd,tau_nd,pos3p_nd,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_radar(mu_nd,t0_1,t_1_nd,ele,posnp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = ele
        self.df = ele_df
        self.rms = rms
        self.method = 'Gibbs/Herrick-Gibbs' 
          
        return self

    def ellipse(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from radar range+angle measurements using Elliptical Orbit Fitting method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test/radar_obs.dat',dtype=str,skiprows=1) 
            >>> # Extract necessary data for IOD from radar range+angle measurements
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> orbele_site = obs_data[:,1:7].astype(float) # Orbital elements(a,ecc,inc,raan,argp,true_anomaly) of the site
            >>> xyz_site = obs_data[:,7:10].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> azalt = obs_data[:,10:12].astype(float) # Azimuth and Altitude angle of space object, [deg]
            >>> r = obs_data[:,12].astype(float) # Slant distance of the space object relative to the site, [km]
            >>> # Load the necessary data to ArcObs
            >>> from orbdtools import ArcObs
            >>> arc_radar = ArcObs({'t':t,'azalt':azalt,'r':r,'xyz_site':xyz_site,'orbele_site':orbele_site}) 
            >>> arc_radar.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_radar.iod(earth)
            >>> arc_iod.ellipse(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of ascending node, [rad] or [deg]
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        """
        if self.mode == 'optical':
            raise Exception('The method of Elliptical Orbit Fitting is not applicable to optical angle-only measurement data for initial orbit determination.')

        posnp_nd = self.posnp_nd
        tof_nd = self.tof_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd

        ele = ellipse_fitting(mu_nd,tof_nd,t_1_nd,posnp_nd,degrees)
 
        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_radar(mu_nd,t0_1,t_1_nd,ele,posnp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = ele
        self.df = ele_df  
        self.rms = rms  
        self.method = 'Elliptical Orbit Fitting'  

        return self

    def circular(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Near-Circular Orbit Hypothesis method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.circular(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of ascending node, [rad] or [deg]
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references: 
            张晓祥,吴连大,熊建宁.空间目标的圆轨道跟踪法[J].天文学报, 2003, 44(4):11.DOI:10.3321/j.issn:0001-5245.2003.04.010.        
        """
        if self.mode == 'radar':
            raise Exception('The method of Near-Circular Orbit Hypothesis is only applicable to angle-only measurement data for initial orbit determination.')

        los3p = self.los3p
        xyz_site3p_nd = self.xyz_site3p_nd

        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        tof_nd = self.tof_nd
        tau_nd = self.tau_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd
        
        ele = solve_circular(mu_nd,tof_nd,tau_nd,los3p,xyz_site3p_nd,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_optical(mu_nd,t0_1,t_1_nd,[ele],losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = ele
        self.df = ele_df   
        self.rms = rms   
        self.method = 'Near-Circular Orbit Hypothesis' 

        return self      

    def gauss(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Gauss method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.gauss(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of ascending node, [rad] or [deg]
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references: 
            Curtis H D. Orbital Mechanics for Engineering Students: Revised 4th edition[M]. Butterworth-Heinemann, 2020. 
        """
        if self.mode == 'radar':
            raise Exception('The method of Gauss is only applicable to angle-only measurement data for initial orbit determination.')

        los3p = self.los3p
        xyz_site3p_nd = self.xyz_site3p_nd

        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        tof_nd = self.tof_nd
        tau_nd = self.tau_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd
        
        eles = gauss_estimate(mu_nd,tof_nd,tau_nd,los3p,xyz_site3p_nd,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_optical(mu_nd,t0_1,t_1_nd,eles,losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = eles
        self.df = ele_df
        self.rms = rms
        self.method = 'Gauss'

        return self   

    def laplace(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Laplace method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.laplace(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of ascending node, [rad] or [deg]
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references:
            Vallado D. Fundamentals of Astrodynamics and Applications(4th)[M], Microcosm Press, 2013. 
            Bate R R, Mueller D D, White J E, et al. Fundamentals of astrodynamics(2nd)[M]. Courier Dover Publications, 2020.          
        """
        if self.mode == 'radar':
            raise Exception('The method of Laplace is only applicable to angle-only measurement data for initial orbit determination.')

        los3p = self.los3p
        xyz_site3p_nd = self.xyz_site3p_nd

        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        tof_nd = self.tof_nd
        tau_nd = self.tau_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd
        
        eles = laplace_estimate(mu_nd,tof_nd,tau_nd,los3p,xyz_site3p_nd,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_optical(mu_nd,t0_1,t_1_nd,eles,losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = eles
        self.df = ele_df
        self.rms = rms
        self.method = 'Laplace'

        return self  

    def multilaplace(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Intermediate epoch from optical angle-only measurement data using Multiple-points Laplace method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.multilaplace(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of asc
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references:
            Bate R R, Mueller D D, White J E, et al. Fundamentals of astrodynamics(2nd)[M]. Courier Dover Publications, 2020.
        """
        if self.mode == 'radar':
            raise Exception('The method of Multiple-points Laplace is only applicable to angle-only measurement data for initial orbit determination.')

        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        t0_2 = self.t0_2
        t_2_nd = self.t_2_nd
        mu_nd = self.mu_nd
        
        eles = multilaplace_estimate(mu_nd,t_2_nd,xyz_sitenp_nd,losnp,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_optical(mu_nd,t0_2,t_2_nd,eles,losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = eles
        self.df = ele_df
        self.rms = rms
        self.method = 'Multiple-points Laplace'
            
        return self 

    def doubleR(self,degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Double-R method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.doubleR(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of asc
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references: 
            Vallado D. Fundamentals of Astrodynamics and Applications(4th)[M], Microcosm Press, 2013. 
        """
        if self.mode == 'radar':
            raise Exception('The method of Double-R is only applicable to angle-only measurement data for initial orbit determination.')

        los3p = self.los3p
        xyz_site3p_nd = self.xyz_site3p_nd

        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        tof_nd = self.tof_nd
        tau_nd = self.tau_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd
        
        ele = doubleR_estimate(mu_nd,tof_nd,tau_nd,los3p,xyz_site3p_nd,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_optical(mu_nd,t0_1,t_1_nd,[ele],losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = ele
        self.df = ele_df
        self.rms = rms
        self.method = 'Double-R'

        return self         

    def gooding(self,tm=1,M=0,method='universal',degrees=True,ellipse_only=True,rms_tol=2e-4):
        """
        Estimate the classical orbital elements at Median epoch from optical angle-only measurement data using Gooding method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.gooding(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            tm -> [int,optional,default=1] Lambert transfer mode. If tm = 1, then short way transfer mode; else if tm = -1, long way transfer mode.  
            M -> [int,optional,default=0] Number of full revolutions in transfer
            method -> [str,optional,default='universal'] Method for solving the Lambert's problem. Available options include 'universal' and 'izzo'
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of asc
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        References:
            Gooding R H. A new procedure for the solution of the classical problem of minimal orbit determination from three lines of sight[J]. Celestial Mechanics and Dynamical Astronomy, 1996, 66: 387-423.
            Izzo D. Revisiting Lambert’s problem[J]. Celestial Mechanics and Dynamical Astronomy, 2015, 121: 1-15.
            Curtis H D. Orbital Mechanics for Engineering Students: Revised 4th edition[M]. Butterworth-Heinemann, 2020. 
        """
        if self.mode == 'radar':
            raise Exception('The method of Gooding is only applicable to angle-only measurement data for initial orbit determination.')

        los3p = self.los3p
        xyz_site3p_nd = self.xyz_site3p_nd

        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        tof_nd = self.tof_nd
        tau_nd = self.tau_nd
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd
        
        ele = gooding_estimate(mu_nd,tof_nd,tau_nd,los3p,xyz_site3p_nd,tm,M,method,degrees)

        # check the orbital elements by rms
        ele_df,rms = to_ele_dict_optical(mu_nd,t0_1,t_1_nd,[ele],losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = ele
        self.df = ele_df
        self.rms = rms
        self.method = 'Gooding'

        return self 

    def fg_series(self,degrees=True,ellipse_only=True,rms_tol=2e-4,improved=True):
        """
        Estimate the classical orbital elements at Intermediate epoch from optical angle-only measurements or radar range+angle measurements using FG-Series method.

        Usage:
            >>> # Load the observation file
            >>> import numpy as np
            >>> obs_data = np.loadtxt('test5/T25872_KUN2_2.dat',dtype=str,skiprows=1) # Load the observation file
            >>> obs_data = obs_data[::10]
            >>> # extract the necessary data
            >>> t = obs_data[:,0] # Obsevation time in UTC
            >>> radec = obs_data[:,1:3].astype(float) # Ra and Dec of space object, [hour,deg]
            >>> xyz_site = obs_data[:,3:6].astype(float) # Cartesian coordinates of the site in GCRF, [km]
            >>> radec[:,0] *= 15 # Convert hours to degrees
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site}) # Load the necessary data
            >>> arc_optical.lowess_smooth() # Eliminate outliers
            >>> # Set the Earth as the central body of attraction
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
            >>> arc_iod = arc_optical.iod(earth)
            >>> arc_iod.fg_series(ellipse_only=False)
            >>> print(arc_iod.df)
        Inputs:
            degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
            ellipse_only -> [bool,optional,default=True] Switch for filtering out elliptical orbits with semi-major axis greater than the equatorial radius of the central body and O-C less than the preset threshold. 
            If True, only elliptical orbits are valid, otherwise all orbits including parabolic and hyperbolic orbits are also valid.
            rms_tol -> [float,optional,default=2e-4] Tolerance of RMS of O-C. If rms > rms_tol, the orbital elements is considered invalid, othwewise valid.
        Outputs:
            ele_df -> [Dataframe of Pandas] Dataframe of possible classical orbital elements with keys and values as follows
                epoch -> [str] Epoch of orbital elements in UTC
                a -> [float] Semi-major axis in [L_nd]
                ecc -> [float] Eccentricity
                inc -> [float] Inclination, [rad] or [deg]
                raan -> [float] Longitude of asc
                argp -> [float] Argument of perigee, [rad] or [deg]
                nu -> [float] True anomaly, [rad] or [deg]
                M -> [float] Mean anomaly, [rad] or [deg]
                h -> [float] Modulus of angular momentum in [L_nd]^2/[T_nd]
                status -> [str] Status of IOD determined by RMS of O-C: 'success' or 'failed'. 
                The status of IOD are only for reference, meaning that although it shows success, the orbital elements may converge to a trivial or even wrong orbit.
        references: 
            李光宇.天体测量和天体力学基础[M].科学出版社,2015.   
            刘林.卫星轨道力学算法[M].南京大学出版社,2019.
            Bate R R, Mueller D D, White J E, et al. Fundamentals of astrodynamics(2nd)[M]. Courier Dover Publications, 2020.     
        Note:
            For radar measurement type point data (3D Cartesian coordinates of space targets), 
            using the IOD results of the Gibbs/Herrick-Gibbs method as the initial value of the FG-Series method 
            cannot effectively improve the results of directly using the FG-Series method.
        """
        t0_1 = self.t0_1
        t_1_nd = self.t_1_nd
        mu_nd = self.mu_nd

        if self.mode == 'radar':
            posnp_nd = self.posnp_nd
            ele = fg_series_radar(mu_nd,t_1_nd,posnp_nd,degrees)
            # check the orbital elements by rms
            ele_df,rms = to_ele_dict_radar(mu_nd,t0_1,t_1_nd,ele,posnp_nd,degrees,rms_tol)
            
        elif self.mode == 'optical':    
            losnp = self.losnp
            xyz_sitenp_nd = self.xyz_sitenp_nd

            if improved:
                self.circular(degrees,ellipse_only,rms_tol)
                coe = self.coe
                rv0 = coe2rv(coe,mu_nd,degrees)
            else:
                rv0 = None    

            ele = fg_series_optical(mu_nd,t_1_nd,xyz_sitenp_nd,losnp,degrees,rv0)
            # check the orbital elements by rms
            ele_df,rms = to_ele_dict_optical(mu_nd,t0_1,t_1_nd,[ele],losnp,xyz_sitenp_nd,degrees,rms_tol)

        if ellipse_only:
            perigee = ele_df['a']*(1-ele_df['ecc'])
            valid_flag = (perigee > self.Re_nd) & (ele_df['status'] == 'success')
            ele_df = ele_df[valid_flag]

        self.coe = ele
        self.df = ele_df
        self.rms = rms
        self.method = 'FG-Series'

        return self  

    '''   
    # Considering that the method of Reference Vector and the method from karimi cannot successfully determine the orbit, the code related to this method is temporarily commented out.
    def ref_vec(self,degrees=True):
        """
        Reference Vector method

        reference:
            贾沛璋,吴连大.初轨计算的参考矢量法[J].天文学报, 1997, 38(4):6.DOI:CNKI:SUN:TWXB.0.1997-04-001.
        """

        if self.mode == 'radar':
            raise Exception('The method of Double-R is only applicable to angular measurement data for initial orbit determination.')
            
        los3p = self.los3p
        xyz_site3p_nd = self.xyz_site3p_nd    

        ta = self._ta_lowess 
        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd

        tof = self.tof
        alpha_vec = self._alpha_vec
        delta_vec = self._delta_vec
        
        ele = ref_vec_estimate(xyz_site3p_nd,los3p,xyz_sitenp_nd,losnp,tof,ta,alpha_vec,delta_vec,degrees)

        # check the orbital elements
        tn = (ta[-1] - ta[0])/2 + ta[0] # Epoch at which the orbit elements is estimated

        ele_dict = to_ele_dict_optical(tn,ta,ele,losnp,xyz_sitenp_nd,degrees)

        return ele_dict   

    def karimi(self,degrees=True):
        """
        reference:
            Karimi R R, Mortari D. Initial orbit determination using multiple observations[J]. Celestial Mechanics and Dynamical Astronomy, 2011, 109: 167-180.
        """

        if self.mode == 'radar':
            raise Exception('The method of Double-R is only applicable to angular measurement data for initial orbit determination.')

        ta = self._ta_lowess 
        losnp = self.losnp
        xyz_sitenp_nd = self.xyz_sitenp_nd
        
        tn,ele = karimi_estimate(ta,xyz_sitenp_nd,losnp,degrees)

        ele_dict = to_ele_dict_optical(tn,ta,ele,losnp,xyz_sitenp_nd,degrees)

        return ele_dict    

    '''                         