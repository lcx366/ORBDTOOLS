import numpy as np
import pandas as pd
from astropy.time import Time
from skyfield.api import EarthSatellite
from sgp4.api import Satrec,SatrecArray
from copy import deepcopy
from os import makedirs,mkdir,path,remove
from spacetrack import SpaceTrackClient
from pathlib import Path
from glob import glob
from datetime import datetime
from colorama import Fore
from time import sleep

from ..utils import data_prepare
from ..arcmatch.parse_tle import load_tle_file
from ..transform.frame_trans import gcrf_teme_mat

class TLE(object):
    """
    class of Two Line Elements
        Attributes:
            - df -> TLE database in form of pandas dataframe, which includes NORAD IDs, mean orbital elements(a, ecc, inc, raan, argp, M), h, bstar, and epoch.
            - range_epoch -> Epoch coverage of the TLE file in format of [epoch_min, epoch_max]
            - h_uec -> Direction vector of angular momentum
            - adot -> Rate of the semi-major axis in [L_nd/min] at the refenence epoch, where L_nd = 6378.135 km for Reference Earth Model - WGS72
            - edot -> Rate of the orbital eccentricity at the refenence epoch
            - raandot -> Rate of the ascending node right ascension in [rad/min] at the refenence epoch
            - argpdot -> Rate of the peripheral angular distance in [rad/min] at the refenence epoch
            - n -> Mean motion in [rad/min] at the refenence epoch
            - ndot -> The first derivative of the mean motion in [rad/min^2] at the refenence epoch
            - nddot -> The second derivative of the mean motion in [rad/min^3] at the refenence epoch
            - _sats_Satrec -> list of space objects in class of Satrec
            - _sats_EarthSatellite -> list of space objects in class of EarthSatellite
            - _statistic -> Statistics of the TLE database
        Methods:
            - retrieve -> Retrieve space objects from the TLE database according to the NORAD IDs.
            - atEpoch -> Calculate the mean orbital elements at a certain epoch.
            - predict -> Calculate the cartesian coordinates of space objects in GCRF(Geocentric Celetial Reference Frame) over a period of time.    
            - download -> Download TLE files from [SPACETRACK](https://www.space-track.org)
    """
    def __init__(self,sats_Satrec,sats_EarthSatellite,info):
        """
        Initialize an instance of class TLE.
        """

        df = info['df']
        df_statistic = df[['a','ecc','inc','raan','argp','h','bstar','mjd']].describe().loc['mean':]
        mjd_statistic = df_statistic['mjd'].loc[['mean','min','25%','50%','75%','max']]
        df_statistic['mjd'].loc[['mean','min','25%','50%','75%','max']] = Time(mjd_statistic, format='mjd').isot
        mjd_statistic = df_statistic['mjd']
        range_epoch = mjd_statistic.loc[['min','max']].to_list()
        df_statistic.rename(columns={"mjd": "epoch"},inplace=True)

        for key in info.keys():
            setattr(self, key, info[key])

        self.range_epoch = range_epoch
        self._statistic = df_statistic
        self._sats_Satrec = sats_Satrec
        self._sats_EarthSatellite = sats_EarthSatellite

    def __repr__(self):
        """
        Returns a more information-rich string representation of the TLE object.
        """
        counts = len(self.adot)
        statistic = self._statistic
        return '<TLE object: Counts = {:d} Statistic = \n{}>'.format(counts,self._statistic)    

    def from_file(fn,t=None,out_days=7):
        """
        Load a TLE file, and return an instance of class TLE.

        Usage:
            >>> from orbdtools import TLE
            >>> tle_file = 'test/tle_20220524.txt'
            >>> epoch_obs = '2022-05-24T08:38:34.000Z' # Approximate epoch of the observed arc, which is optional
            >>> tle = TLE.from_file(tle_file,epoch_obs) # Only TLEs within a week before and after the observation epoch are considered
            >>> # tle = TLE.from_file(tle_file) # All TLEs are considered
        Inputs:
            fn -> [str/list of str] Filename of the TLE file, such as 'test/tle_20220524.txt'; or list of str, such as ['line1','line2']
            t -> [str,optional,default=None] Approximate epoch of the observed arc
            out_days -> [float,optional,default=7] Valid time interval of TLE. Only TLEs within out_days before and after the observation epoch are considered
        Outputs:
            tle -> Instance of class TLE
        """
        ts = data_prepare.ts

        # load and parse the TLE data
        fn_type = type(fn)
        if fn_type is str:
            sats_Satrec,sats_EarthSatellite = load_tle_file(fn)
        elif fn_type is list:
            line1,line2 = fn
            satrec = Satrec.twoline2rv(line1, line2)
            sats_Satrec = np.array([satrec])
            sats_EarthSatellite = np.array([EarthSatellite.from_satrec(satrec, ts)])
            
        data,h_uec = [],[]
        adot,edot,raandot,argpdot = [],[],[],[]
        n,ndot,nddot = [],[],[]
        unvalid_index = []

        if t is not None: 
            ta = Time(t)
            tsf = ts.from_astropy(ta)

        j = 0
        for sat in sats_EarthSatellite:
            epoch_jd = sat.epoch # in time scale of TT
            if t is not None:
                if np.abs(epoch_jd - tsf) > out_days: 
                    unvalid_index.append(j)
                    continue
            noradid = sat.model.satnum
            a = sat.model.a # normalized semi-major axis
            ecc = sat.model.ecco
            inc_rad,raan_rad,argp_rad,M_rad = sat.model.inclo,sat.model.nodeo,sat.model.argpo,sat.model.mo # in rad
            inc_deg,raan_deg,argp_deg,M_deg = np.rad2deg(inc_rad),np.rad2deg(raan_rad),np.rad2deg(argp_rad),np.rad2deg(M_rad) # in deg
            h = np.sqrt(a*(1-ecc**2))
            bstar = sat.model.bstar
            epoch = epoch_jd.utc_iso(places=3)
            mjd = Time(epoch).mjd
            data.append([noradid,a,ecc,inc_deg,raan_deg,argp_deg,M_deg,h,bstar,epoch,mjd])
            h_uec.append([np.sin(raan_rad)*np.sin(inc_rad),-np.cos(raan_rad)*np.sin(inc_rad),np.cos(inc_rad)]) # unit vector of Angular momentum

            # approximate rate
            raandot.append(sat.model.nodedot)
            argpdot.append(sat.model.argpdot) # in rad/min
            # note that sat.model.ndot = half value of the first derivative of the mean motion
            # sat.model.nddot = Second derivative of the mean motion divided by 6
            n.append(sat.model.no)
            ndot.append(sat.model.ndot*2)
            nddot.append(sat.model.nddot*6) # in rad/min, rad/min**2, rad/min**3
            adot.append(-2/3*a*sat.model.ndot/sat.model.no) # in L_nd/min
            edot.append(-2/3*(1-ecc)*sat.model.ndot/sat.model.no) # in 1/min  

            j += 1

        h_uec = np.array(h_uec)
        adot,edot,raandot,argpdot = np.array(adot),np.array(edot),np.array(raandot),np.array(argpdot)
        n,ndot,nddot = np.array(n),np.array(ndot),np.array(nddot)   
    
        # Statistics for TLE database
        df = pd.DataFrame(data, columns=['noradid','a','ecc','inc','raan','argp','M','h','bstar','epoch','mjd'])

        info = {'df':df,'h_uec':h_uec,'adot':adot,'edot':edot,'raandot':raandot,'argpdot':argpdot,'n':n,'ndot':ndot,'nddot':nddot} 

        sats_Satrec_valid = np.delete(sats_Satrec,unvalid_index) 
        sats_EarthSatellite_valid = np.delete(sats_EarthSatellite,unvalid_index) 

        return TLE(sats_Satrec_valid,sats_EarthSatellite_valid,info)

    def retrieve(self,sats_id=None):
        """
        Retrieve space objects from the TLE database according to their NORAD IDs.

        Usage:
            >>> sats_id = [47,58,52139,52140,52150] # NORAD IDs of space objects
            >>> tle_retrieve = tle.retrieve(sats_id)
        Inputs:
            sats_id -> [list of int,optinal,default=None] NORAD IDs of space objects. If None, all space objects will participate in the calculation.  
        Outputs:
            tle_retrieve -> Instance of class TLE
        """
        if sats_id is None:
            in_flag = np.ones_like(self.n).astype(bool)
        else:
            in_flag = self.df['noradid'].isin(sats_id)

        sats_Satrec = np.array(self._sats_Satrec)[in_flag]
        sats_EarthSatellite = np.array(self._sats_EarthSatellite)[in_flag]   
        
        info = deepcopy(self.__dict__)

        for key in info.keys():
            info[key] = info[key][in_flag]

        return TLE(sats_Satrec,sats_EarthSatellite,info)

    def atEpoch(self,epoch,sats_id=None):
        """
        Calculate the mean orbital elements at a certain epoch.

        Usage:
            >>> epoch_obs = '2022-05-24T08:38:34.000Z'
            >>> tle_epoch = tle.atEpoch(epoch_obs)
            >>> # sats_id = [10,47,58,52139]
            >>> # tle_epoch = tle.atEpoch(epoch_obs,sats_id)
        Inputs:
            epoch -> [str or Astropy Time object] Epoch at which to calculate the mean orbital elements
            sats_id -> [list of int,optinal,default=None] NORAD IDs of space objects. If None, all space objects will participate in the calculation.  
        Outputs:
            tle_epoch -> Instance of class TLE with updated mean orbital elements
        """

        if sats_id is None:
            in_flag = np.ones_like(self.n).astype(bool)
        else:
            in_flag = self.df['noradid'].isin(sats_id)

        sats_Satrec = np.array(self._sats_Satrec)[in_flag]
        sats_EarthSatellite = np.array(self._sats_EarthSatellite)[in_flag]   
        
        info = deepcopy(self.__dict__)

        for key in info.keys():
            info[key] = info[key][in_flag]

        df = info['df']
        h_uec = info['h_uec']

        ts = data_prepare.ts
        ta = Time([epoch]*len(sats_Satrec))
        tsf = ts.from_astropy(ta) 

        j = 0
        for sat,t_i in zip(sats_EarthSatellite,tsf):
            sat.at(t_i)
            inc_rad = sat.model.im
            raan_rad = sat.model.Om

            a = sat.model.am
            ecc = sat.model.em
            inc = np.rad2deg(inc_rad)
            raan = np.rad2deg(raan_rad) % 360
            argp = np.rad2deg(sat.model.om) % 360
            M = np.rad2deg(sat.model.mm) % 360
            h = np.sqrt(sat.model.am*(1-sat.model.em**2))

            df.iloc[j,1:-3] = [a,ecc,inc,raan,argp,M,h]

            h_uec[j] = [np.sin(raan_rad)*np.sin(inc_rad),-np.cos(raan_rad)*np.sin(inc_rad),np.cos(inc_rad)] # unit vector of Angular momentum

            j += 1   

        df['epoch'] = ta.isot
        df['epoch'] += 'Z'
        df['mjd'] = ta.mjd

        return TLE(sats_Satrec,sats_EarthSatellite,info)    
            
    def predict(self,t,sats_id=None):
        """
        Calculate the cartesian coordinates of space objects in GCRF(Geocentric Celetial Reference Frame) over a period of time.

        Usage:
            >>> t_list = ['2022-05-23T16:22:49.408Z', '2022-05-23T18:48:34.488Z',
                          '2022-05-23T18:09:35.640Z', '2022-05-23T20:54:20.228Z',
                          '2022-05-23T21:10:37.703Z', '2022-05-23T17:48:24.865Z']
            >>> xyz_gcrf = tle.predict(t_list)
        Inputs:
            t -> [str or list of str] time to perform a orbital propagation
            sats_id -> [list of int,optinal,default=None] NORAD IDs of space objects  
        Outputs:
            xyz_gcrf -> [3D array] Cartesian coordinates of space objects in GCRF    
        """
        ta = Time(t)

        if sats_id is None:
            sats = self._sats_Satrec
        else:
            in_flag = self.df['noradid'].isin(sats_id)
            sats = np.array(self._sats_Satrec)[in_flag]

        sats_array = SatrecArray(sats)  

        if ta.ndim == 0:
            e, xyz_teme, vxyz_teme = sats_array.sgp4(np.array([ta.jd1]),np.array([ta.jd2]))
        else:    
            e, xyz_teme, vxyz_teme = sats_array.sgp4(ta.jd1,ta.jd2)

        gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta)
        xyz_gcrf = (teme2gcrf_mat @ xyz_teme.transpose(1,2,0)).transpose(2,0,1)  

        return xyz_gcrf   

    def download(noradids,mode='keep',dir_TLE='TLE/'):
        """
        Download TLE files from [SPACETRACK](https://www.space-track.org)

        Usage: 
            >>> tlefile = TLE.download(noradids)
            >>> tlefile = TLE.download(noradids,'clear')
            >>> tlefile = TLE.download('satno.txt')
        Inputs:
            noradids -> [str, int, list of str/int] NORADID of space targets. It can be a single NORADID, list of NORADID, or a file containing a set of NORADIDs.
            A sample of the file is as follows:
            #satno
            12345
            23469
            45678
            mode -> [str,optional,default='keep'] If 'keep', files stored in TLE directory will be intact during the downloading; 
            if 'clear', files stored in TLE directory will be removed ahead of downloading.
            dir_TLE -> [str,optional,default='TLE/'] Directory of the TLE file to store
        Outputs: 
            tlefile  -> [str] Path of the TLE file downloaded
        """
        # Check whether a list is empty or not
        if not noradids: raise Exception('noradids is empty.')

        if type(noradids) is list:
            if type(noradids[0]) is int: noradids = [str(i) for i in noradids]    
        else:
            noradids = str(noradids)
            if '.' in noradids: # noradids as a file
                noradids = list(set(np.loadtxt(noradids,dtype=str)))
            else:
                noradids = [noradids]    
    
        # Set the maximum of requested URL's length with a single access 
        # The setup prevents exceeding the capacity limit of the server
        n = 500
        noradids_parts = [noradids[i:i + n] for i in range(0, len(noradids), n)]  
        part_num = len(noradids_parts)    
    
        # username and password for Space-Track
        home = str(Path.home())
        direc = home + '/src/spacetrack-data/'
        loginfile = direc + 'spacetrack-login'

        if not path.exists(direc): makedirs(direc)
        if not path.exists(loginfile):
            username = input('Please input the username for Space-Track(which can be created at https://www.space-track.org/auth/login): ')
            password = input('Please input the password for Space-Track: ')
            outfile = open(loginfile,'w')
            for element in [username,password]:
                outfile.write('{:s}\n'.format(element))
            outfile.close()
        else:
            infile = open(loginfile,'r')
            username = infile.readline().strip()
            password = infile.readline().strip()
            infile.close()
    
        # save TLE data to files  
        fileList_TLE = glob(dir_TLE+'*')
        if path.exists(dir_TLE):
            if mode == 'clear':
                for file in fileList_TLE:
                    remove(file)
        else:
            makedirs(dir_TLE) 

        valid_ids,j = [],1
        date_str = datetime.utcnow().strftime("%Y%m%d")
        tlefile = dir_TLE + 'tle_{:s}.txt'.format(date_str)
        file_tle = open(tlefile,'w')  

        st = SpaceTrackClient(username, password)
        for part in noradids_parts:
            desc = 'Downloading TLE data: Part {:s}{:2d}{:s} of {:2d}'.format(Fore.BLUE,j,Fore.RESET,part_num)
            print(desc,end='\r')

            try:
                lines_tle = st.tle_latest(norad_cat_id=part,ordinal=1,iter_lines=True,format='tle') 
            except:       
                remove(loginfile)
                raise ConnectionError("401 Unauthorized: username or password entered incorrectly!") 
                    
            for line in lines_tle:
                words = line.split()
                if words[0] == '2': valid_ids.append(words[1].lstrip('0'))
                file_tle.write(line+'\n')
            sleep(j+5) 
            j += 1   
        file_tle.close()

        missed_ids = list(set(noradids)-set(valid_ids))
        if missed_ids: 
            missed_ids_filename = dir_TLE + 'missed_ids_{:s}.txt'.format(date_str)
            desc = '{:s}Note: space targets with unavailable TLE are stored in {:s}.{:s} '.format(Fore.RED,missed_ids_filename,Fore.RESET)
            print(desc) 
            np.savetxt(missed_ids_filename,missed_ids,fmt='%s')

        return tlefile