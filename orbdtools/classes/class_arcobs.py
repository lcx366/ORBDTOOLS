from astropy.time import Time

from ..arcmatch.optical import arcsat_match as arcsat_match_optical
from ..arcmatch.radar import arcsat_match as arcsat_match_radar
from ..utils.preprocessing import lowess_smooth_optical,lowess_smooth_radar

class ArcObs(object):
    """
    class of Observation Arc
        Attributes:
            - mode -> [str] Type of the observation arc. 'optical' or 'radar' are avaliable.
            - t -> [array of str] Time sequence from the observation arc
            - radec -> [2D array] Ra and Dec of space objects from the observation arc, [deg]
            - xyz_site -> [3D array] Cartesian coordinates of space objects in GCRF from the observation arc, [km]
            - r -> [array of float] Slant distance of the space object relative to the site, [km]
            - _ta_lowess -> [array like] Astropy Time sequence with outliers removed using LOWESS 
            - _radec_lowess -> [2D array] Ra and Dec with outliers removed using LOWESS 
            - _xyz_site_lowess -> [3D array] Cartesian coordinates with outliers removed using LOWESS 
            - _r_lowess -> [array of float] Slant distance with outliers removed using LOWESS 
            - flag_lowess -> [array of bool] Flag of data points. If False, the data point is an outlier.
        Methods:    
            - lowess_smooth -> 
            - arc_match
    """

    def __init__(self,info):
        """
        Create an instance of class ArcObs

        Usage:
            >>> from orbdtools import ArcObs
            >>> arc_optical = ArcObs({'t':t,'radec':radec,'xyz_site':xyz_site})
            >>> # arc_radar = ArcObs({'t':t,'radec':radec,'r':r,'xyz_site':xyz_site})
        Inputs:
            info -> [dict] Necessary data for observation arcs in form of dictionary    
        Outputs:
            arc_optical -> An instance of class ArcObs 
        """
        if 'r' in info.keys():
            info['mode'] = 'radar'
        else:
            info['mode'] = 'optical'

        info['_ta_lowess'] = Time(info['t'])
        info['_radec_lowess'] = info['radec']
        info['_xyz_site_lowess'] = info['xyz_site'] 
        if 'r' in info.keys(): info['_r_lowess'] = info['r']

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
        elif self.mode == 'radar':
            flag = lowess_smooth_radar(ta,self.radec,self.r,frac)
            self._r_lowess = info['_r_lowess'] = self.r[flag]

        self._ta_lowess = info['_ta_lowess'] = ta[flag]
        self._radec_lowess = info['_radec_lowess'] = self.radec[flag]
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
        radec = self._radec_lowess
        xyz_site = self._xyz_site_lowess

        # For case of optical angle measurement data
        if self.mode == 'optical':
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
            r = self._r_lowess
            if threshold_dict is None:
                threshold_pre = [10,10] # [deg,km]
                threshold_deep = [2000,5] # [arcsec,km]
                threshold_slope = [400,2] # [arcsec/s,km/s]
            else:
                threshold_pre = threshold_dict['threshold_pre']
                threshold_deep = threshold_dict['threshold_deep']
                threshold_slope = threshold_dict['threshold_slope'] 
            code_match,satnum,disp_match = arcsat_match_radar(tle,ta,xyz_site,radec,r,threshold_pre,threshold_deep,threshold_slope)

        self.code_match = info['code_match'] = code_match
        self.satnum = info['satnum'] = satnum
        self.disp_match = info['disp_match'] = disp_match  

        return None