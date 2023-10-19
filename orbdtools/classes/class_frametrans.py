from ..transform import frame_trans                                              

class MatrixTrans(object):
    
    def __init__(self,info):

        for key in info.keys():
            setattr(self, key, info[key])   

    def __repr__(self):
        """
        Returns a more information-rich string representation of the MatrixTrans object.
        """
        frame1,frame2 = self._transname.split('_')
        return "<MatrixTrans object: '{:s}' ⇌ '{:s}'>".format(frame1.upper(),frame2.upper())      

class FrameTrans(object):
    """
    Class of reference frame transformation.

        Methods:    
            - lrf_topo_mat -> Calculate the transformation matrix between the Launch Reference Frame and the Topocentric NEU(North East Up) Reference Frame.
            - topo_itrf_mat -> Calculate the transformation matrix between the Topocentric NEU(North East Up) and the ITRF(International Terrestrial Reference Frame).
            - gcrf_itrf_mat -> Calculate the transformation matrix between the GCRF(Geocentric Celestial Reference Frame) and the ITRF(International Terrestrial Reference Frame).
            - gcrf_topo_mat -> Calculate the transformation matrix between the GCRF(Geocentric Celestial Reference Frame) and the Topocentric NEU(North East Up) Reference Frame.
            - gcrf_teme_mat -> Calculate the transformation matrix between the GCRF(Geocentric Celestial Reference Frame) and the TEME(True Equator, Mean Equinox) Reference Frame.
            - meme_topo_mat -> Calculate the transformation matrix between the MEME(Mean Equator, Mean Equinox) Reference Frame and the Topocentric NEU(North East Up) Reference Frame.
            - lrf_gcrf_mat -> Calculate the transformation matrix between the Launch Reference Frame and GCRF(Geocentric Celestial Reference Frame).
            - lirf_lrf_mat -> Calculate the transformation matrix between the Launch Inertial Reference Frame and the Launch Reference Frame.
            - lrf_meme_mat -> Calculate the transformation matrix between the Launch Reference Frame and the MEME(Mean Equator, Mean Equinox) Reference Frame.
            - lrf_teme_mat -> Calculate the transformation matrix between the Launch Reference Frame and the TEME(True Equator, Mean Equinox) Reference Frame.
            - ECI_PQW_mat -> Calculate the transformation matrix between the Earth Centred Inertial (ECI) reference frame and the PQW(perifocal) reference frame.
            - ECI_RSW_mat -> Calculate the transformation matrix between the Earth Centred Inertial (ECI) reference frame and the RSW(x:Radial,y:Along-track,z:Cross-track) reference frame.
            - ECI_NTW_mat -> Calculate the transformation matrix between the Earth Centred Inertial (ECI) reference frame and the NTW(x:Normal,y:Tangent,z:Cross-track) reference frame.
            - ECI_RADAR_mat -> Calculate the transformation matrix between the Earth Centred Inertial (ECI) reference frame and the RADAR(x:Along-track,y:Cross-track,z:Radial) reference frame.
            - RSW_BF_mat -> Calculate the transformation matrix between the RSW(x:Radial,y:Along-track,z:Cross-track) reference frame and the Body-Fixed reference frame.
            - BF_DF_mat -> Calculate the transformation matrix between the Body-Fixed(BF) reference frame and the Device-Fixed(DF) reference frame.
            - ECI_DF_mat -> Calculate the transformation matrix between the ECI reference frame and the Device-Fixed reference frame.
    """
    def lrf_topo_mat(alpha,degrees=True): 
        """
        Rotation matrix between the Launch Reference Frame and the Topocentric NEU(North East Up) Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> alpha = 30 # Azimuth of the launch direction in [deg]
            >>> #alpha = [30,40]
            >>> matrix_trans = FrameTrans.lrf_topo_mat(alpha)
            >>> lrf2topo_mat = matrix_trans.lrf2topo_mat
            >>> topo2lrf_mat = matrix_trans.topo2lrf_mat
        Inputs:
            alpha -> [float,list of float] Azimuth of the launch direction in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, the azimuth is assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        lrf2topo_mat,topo2lrf_mat = frame_trans.lrf_topo_mat(alpha,degrees)
        info = {'lrf2topo_mat':lrf2topo_mat,'topo2lrf_mat':topo2lrf_mat,'_transname':'lrf_topo'}
        return MatrixTrans(info)

    def topo_itrf_mat(lon,lat,degrees=True): 
        """
        Rotation matrix between the Topocentric NEU(North East Up) and the ITRF(International Terrestrial Reference Frame).

        Usage:
            >>> from orbdtools import FrameTrans
            >>> lon,lat = 102.5,25.2 # longitude and latitude in [deg]
            >>> #lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> matrix_trans = FrameTrans.topo_itrf_mat(lon,lat) 
            >>> topo2itrf_mat = matrix_trans.topo2itrf_mat
            >>> itrf2topo_mat = matrix_trans.itrf2topo_mat
            >>> print(topo2itrf_mat)
            >>> print(itrf2topo_mat)
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, the longitude and latitude are assumed to be in degrees
        Outputs: 
            matrix_trans -> Instance of class MatrixTrans
        """
        topo2itrf_mat,itrf2topo_mat = frame_trans.topo_itrf_mat(lon,lat,degrees)
        info = {'topo2itrf_mat':topo2itrf_mat,'itrf2topo_mat':itrf2topo_mat,'_transname':'topo_itrf'}
        return MatrixTrans(info) 

    def gcrf_itrf_mat(ta):
        """
        Rotation matrix between the GCRF(Geocentric Celestial Reference Frame) and the ITRF(International Terrestrial Reference Frame).

        Usage:
            >>> from astropy.time import Time
            >>> ta = Time('2022-01-02T03:04:05.678') # UTC
            >>> # ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> from orbdtools import FrameTrans
            >>> matrix_trans = FrameTrans.gcrf_itrf_mat(ta)
            >>> gcrf2itrf_mat = matrix_trans.gcrf2itrf_mat
            >>> itrf2gcrf_mat = matrix_trans.itrf2gcrf_mat
            >>> print(gcrf2itrf_mat)
            >>> print(itrf2gcrf_mat)
        Inputs:
            ta -> [array-like of Astropy Time object] Time to make transformation
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        gcrf2itrf_mat,itrf2gcrf_mat = frame_trans.gcrf_itrf_mat(ta)
        info = {'gcrf2itrf_mat':gcrf2itrf_mat,'itrf2gcrf_mat':itrf2gcrf_mat,'_transname':'gcrf_itrf'}
        return MatrixTrans(info)       

    def gcrf_topo_mat(lon,lat,ta,degrees=True):
        """
        Rotation matrix between the GCRF(Geocentric Celestial Reference Frame) and the Topocentric NEU(North East Up) Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> lon,lat = 102.5,25.2
            >>> ta = Time(['2022-01-02T03:04:05.678']) # UTC
            >>> #lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> #ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.gcrf_topo_mat(lon,lat,ta)
            >>> gcrf2topo_mat = matrix_trans.gcrf2topo_mat
            >>> topo2gcrf_mat = matrix_trans.topo2gcrf_mat
            >>> print(gcrf2topo_mat)
            >>> print(topo2gcrf_mat)
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg] 
            ta -> [array-like of Astropy Time object] Time to make transformation
            degrees -> [bool,optional,default=True] If True, the longitude and latitude are assumed to be in degrees 
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        gcrf2topo_mat,topo2gcrf_mat = frame_trans.gcrf_topo_mat(lon,lat,ta,degrees)
        info = {'gcrf2topo_mat':gcrf2topo_mat,'topo2gcrf_mat':topo2gcrf_mat,'_transname':'gcrf_topo'}
        return MatrixTrans(info)  

    def gcrf_teme_mat(ta):
        """
        Rotation matrix between the GCRF(Geocentric Celestial Reference Frame) and the TEME(True Equator, Mean Equinox) Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> ta = Time('2022-01-02T03:04:05.678') # UTC
            >>> # ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.gcrf_teme_mat(ta)
            >>> gcrf2teme_mat = matrix_trans.gcrf2teme_mat
            >>> teme2gcrf_mat = matrix_trans.teme2gcrf_mat 
            >>> print(gcrf2teme_mat)
            >>> print(teme2gcrf_mat)
        Inputs:
            ta -> [array-like of Astropy Time object] Time to make transformation
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        gcrf2teme_mat,teme2gcrf_mat = frame_trans.gcrf_teme_mat(ta)
        info = {'gcrf2teme_mat':gcrf2teme_mat,'teme2gcrf_mat':teme2gcrf_mat,'_transname':'gcrf_teme'}
        return MatrixTrans(info)        

    def meme_topo_mat(lon,lat,ta,degrees=True): 
        """
        Rotation matrix between the MEME(Mean Equator, Mean Equinox) Reference Frame and the Topocentric NEU(North East Up) Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> lon,lat = 102.5,25.2
            >>> ta = Time('2022-01-02T03:04:05.678') # UTC
            >>> #lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> #ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.meme_topo_mat(lon,lat,ta)
            >>> meme2topo_mat = matrix_trans.meme2topo_mat
            >>> topo2meme_mat = matrix_trans.topo2meme_mat 
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg]
            ta -> [array-like of Astropy Time object] Time to make transformation
            degrees -> [bool,optional,default=True] If True, the longitude and latitude are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        meme2topo_mat,topo2meme_mat = frame_trans.meme_topo_mat(lon,lat,ta,degrees)
        info = {'meme2topo_mat':meme2topo_mat,'topo2meme_mat':topo2meme_mat,'_transname':'meme_topo'}
        return MatrixTrans(info)       

    def lrf_gcrf_mat(lon,lat,alpha,ta,degrees=True):
        """
        Rotation matrix between the Launch Reference Frame and GCRF(Geocentric Celestial Reference Frame).

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> alpha = [20,30,40]
            >>> ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.lrf_gcrf_mat(lon,lat,alpha,ta)
            >>> lrf2gcrf_mat = matrix_trans.lrf2gcrf_mat 
            >>> gcrf2lrf_mat = matrix_trans.gcrf2lrf_mat 
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg]
            alpha -> [float,list of float] Azimuth of the launch direction in [rad] or [deg]
            ta -> [array-like of Astropy Time object] Time to make transformation
            degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        lrf2gcrf_mat,gcrf2lrf_mat = frame_trans.lrf_gcrf_mat(lon,lat,alpha,ta,degrees)
        info = {'lrf2gcrf_mat':lrf2gcrf_mat,'gcrf2lrf_mat':gcrf2lrf_mat,'_transname':'lrf_gcrf'}
        return MatrixTrans(info)  

    def lirf_lrf_mat(lon,lat,alpha,ta0,ta,degrees=True):
        """
        Rotation matrix between the Launch Inertial Reference Frame and the Launch Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> ta0 = Time(['2023-01-02T03:04:05.678','2023-02-02T03:04:05.678']) # UTC
            >>> ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.lirf_lrf_mat(lon,lat,alpha,ta0,ta)
            >>> lirf2lrf_mat = matrix_trans.lirf2lrf_mat
            >>> lrf2lirf_mat = matrix_trans.lrf2lirf_mat 
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg]
            alpha -> [float,list of float] Azimuth of the lrf direction in [rad] or [deg]
            ta0 -> [array-like of Astropy Time object] Epoch at which the inertial reference frame defined
            ta -> [array-like of Astropy Time object] Time to make transformation
            degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        if ta.isscalar ^ ta0.isscalar:
            raise Exception('length of ta0 should be equal to that of ta')
        lirf2lrf_mat,lrf2lirf_mat = frame_trans.lirf_lrf_mat(lon,lat,alpha,ta0,ta,degrees)
        info = {'lirf2lrf_mat':lirf2lrf_mat,'lrf2lirf_mat':lrf2lirf_mat,'_transname':'lirf_lrf'}
        return MatrixTrans(info)   

    def lrf_meme_mat(lon,lat,alpha,ta,degrees=True):
        """
        Rotation matrix between the Launch Reference Frame and the MEME(Mean Equator, Mean Equinox) Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> alpha = [20,30,40]
            >>> ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.lrf_meme_mat(lon,lat,alpha,ta)
            >>> lrf2meme_mat = matrix_trans.lrf2meme_mat 
            >>> meme2lrf_mat = matrix_trans.meme2lrf_mat 
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg]
            alpha -> [float,list of float] Azimuth of the lrf direction in [rad] or [deg]
            ta -> [array-like of Astropy Time object] Time to make transformation
            degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        lrf2meme_mat,meme2lrf_mat = frame_trans.lrf_meme_mat(lon,lat,alpha,ta,degrees)
        info = {'lrf2meme_mat':lrf2meme_mat,'meme2lrf_mat':meme2lrf_mat,'_transname':'lrf_meme'}
        return MatrixTrans(info)  

    def lrf_teme_mat(lon,lat,alpha,ta,degrees=True):
        """
        Rotation matrix between the Launch Reference Frame and the TEME(True Equator, Mean Equinox) Reference Frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> lon,lat = [102.5,102,102.5],[25.2,26,16.5]
            >>> alpha = [20,30,40]
            >>> ta = Time(['2022-01-02T03:04:05.678','2022-02-02T03:04:05.678']) # UTC
            >>> matrix_trans = FrameTrans.lrf_teme_mat(lon,lat,alpha,ta)
            >>> lrf2teme_mat = matrix_trans.lrf2teme_mat 
            >>> teme2lrf_mat = matrix_trans.teme2lrf_mat 
        Inputs:
            lon -> [float,list of float] Longitude of site in [rad] or [deg]
            lat -> [float,list of float] Latitude of site in [rad] or [deg]
            alpha -> [float,list of float] Azimuth of the lrf direction in [rad] or [deg]
            ta -> [array-like of Astropy Time object] Time to make transformation
            degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        lrf2teme_mat,teme2lrf_mat = frame_trans.lrf_teme_mat(lon,lat,alpha,ta,degrees)
        info = {'lrf2teme_mat':lrf2teme_mat,'teme2lrf_mat':teme2lrf_mat,'_transname':'lrf_teme'}
        return MatrixTrans(info)  

    def ECI_PQW_mat(inc,raan,argp,degrees=True):
        """
        Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the PQW(perifocal) reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> inc,raan,argp = 30,40,50 # [i, Ω, ω] in unit of [deg]
            >>> matrix_trans = FrameTrans.ECI_PQW_mat(inc,raan,argp)
            >>> ECI2PQW_mat = matrix_trans.ECI2PQW_mat
            >>> PQW2ECI_mat = matrix_trans.PQW2ECI_mat 
            >>> print(ECI2PQW_mat)
            >>> print(PQW2ECI_mat)
        Inputs:
            inc -> [float, list of float] Orbital inclination (i) in [rad] or [deg]
            raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
            argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        ECI2PQW_mat,PQW2ECI_mat = frame_trans.ECI_PQW_mat(inc,raan,argp,degrees)
        info = {'ECI2PQW_mat':ECI2PQW_mat,'PQW2ECI_mat':PQW2ECI_mat,'_transname':'ECI_PQW'}
        return MatrixTrans(info) 

    def ECI_RSW_mat(inc,raan,argp,nu,degrees=True):
        """
        Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the RSW(x:Radial,y:Along-track,z:Cross-track) reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> inc,raan,argp,nu = 30,40,50,60 # [i, Ω, ω, ν] in unit of [deg]
            >>> matrix_trans = FrameTrans.ECI_RSW_mat(inc,raan,argp,nu)
            >>> ECI2RSW_mat = matrix_trans.ECI2RSW_mat
            >>> RSW2ECI_mat = matrix_trans.RSW2ECI_mat 
            >>> print(ECI2RSW_mat)
            >>> print(RSW2ECI_mat)
        Inputs:
            inc -> [float, list of float] Orbital inclination (i)  in [rad] or [deg]
            raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
            argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
            nu -> [float, list of float] True anomaly (ν) in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        ECI2RSW_mat,RSW2ECI_mat = frame_trans.ECI_RSW_mat(inc,raan,argp,nu,degrees)
        info = {'ECI2RSW_mat':ECI2RSW_mat,'RSW2ECI_mat':RSW2ECI_mat,'_transname':'ECI_RSW'}
        return MatrixTrans(info)     

    def ECI_NTW_mat(ecc,inc,raan,argp,nu,degrees=True):
        """
        Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the NTW(x:Normal,y:Tangent,z:Cross-track) reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> ecc,inc,raan,argp,nu = 0.1,30,40,50,60 # [i, Ω, ω, v] in unit of [deg]
            >>> #ecc,inc,raan,argp,nu = [0.1,0.2],[30,40],[40,50],[50,60],[60,70] # [i, Ω, ω, v] in unit of [deg]
            >>> matrix_trans = FrameTrans.ECI_NTW_mat(ecc,inc,raan,argp,nu)
            >>> ECI2NTW_mat = matrix_trans.ECI2NTW_mat
            >>> NTW2ECI_mat = matrix_trans.NTW2ECI_mat
        Inputs:
            ecc -> [float, list of float] Orbital eccentricity (e)
            inc -> [float, list of float] Orbital inclination (i)  in [rad] or [deg]
            raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
            argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
            nu -> [float, list of float] True anomaly (ν) in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans
        """
        ECI2NTW_mat,NTW2ECI_mat = frame_trans.ECI_NTW_mat(ecc,inc,raan,argp,nu,degrees)
        info = {'ECI2NTW_mat':ECI2NTW_mat,'NTW2ECI_mat':NTW2ECI_mat,'_transname':'ECI_NTW'}
        return MatrixTrans(info) 

    def ECI_RADAR_mat(inc,raan,argp,nu,degrees=True):  
        """
        Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the RADAR(x:Along-track,y:Cross-track,z:Radial) reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> inc,raan,argp,nu = 30,40,50,60 # [i, Ω, ω, v] in unit of [deg]
            >>> #inc,raan,argp,nu = [30,40],[40,50],[50,60],[60,70] # [i, Ω, ω, v] in unit of [deg]
            >>> matrix_trans = FrameTrans.ECI_RADAR_mat(inc,raan,argp,nu)
            >>> ECI2RADAR_mat = matrix_trans.ECI2RADAR_mat
            >>> RADAR2ECI_mat = matrix_trans.RADAR2ECI_mat 
        Inputs:
            inc -> [float, list of float] Orbital inclination (i)  in [rad] or [deg]
            raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
            argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
            nu -> [float, list of float] True anomaly (ν) in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
        Outputs:
            matrix_trans -> Instance of class MatrixTrans 
        """      
        ECI2RADAR_mat,RADAR2ECI_mat = frame_trans.ECI_RADAR_mat(inc,raan,argp,nu,degrees)
        info = {'ECI2RADAR_mat':ECI2RADAR_mat,'RADAR2ECI_mat':RADAR2ECI_mat,'_transname':'ECI_RADAR'}
        return MatrixTrans(info) 

    def RSW_BF_mat(triad,mode,degrees=True): 
        """
        Rotation matrix between the RSW(x:Radial,y:Along-track,z:Cross-track) reference frame and the Body-Fixed reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> triad,mode = [20,30,40],'euler' # Classic 'ZXZ' rotation transform from RSW to BF
            >>> # triad,mode = [20,30,40],'ypr' # Yaw(x)–Pitch(z)–Roll(y) rotation transform from RSW to BF
            >>> # triad,mode = [[0,1,0],[0,0,1],[1,0,0]],'matrix' # Each column of matrix is the base vector of BF in RSW
            >>> # triad,mode = [1,2,3,4],'quaternion' # Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format.
            >>> matrix_trans = FrameTrans.RSW_BF_mat(triad,mode)
            >>> RSW2BF_mat = matrix_trans.RSW2BF_mat
            >>> BF2RSW_mat = matrix_trans.BF2RSW_mat
            >>> print(RSW2BF_mat)
            >>> print(BF2RSW_mat)
            >>> # For multiple dimension
            >>> from orbdtools import FrameTrans
            >>> triad,mode = [[20,30,40],[60,60,70]],'euler' # Classic 'ZXZ' rotation transform from RSW to BF
            >>> #triad,mode = [[20,30,40],[60,60,70]],'ypr' # Yaw(x)–Pitch(z)–Roll(y) rotation transform from RSW to BF
            >>> #triad,mode = [[[0,1,0],[0,0,1],[1,0,0]],[[0,-1,0],[0,0,1],[-1,0,0]]],'matrix' # Each column of matrix is the base vector of BF in RSW
            >>> #triad,mode = [[1,2,3,4],[5,6,7,8]],'quaternion' # Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format.
            >>> matrix_trans = FrameTrans.RSW_BF_mat(triad,mode)
            >>> RSW2BF_mat = matrix_trans.RSW2BF_mat
            >>> BF2RSW_mat = matrix_trans.BF2RSW_mat
        Inputs:
            triad -> [1D/2D array of float] Three angles in [rad] or [deg] for mode of 'euler', 'ypr' and 'quaternion'
                  -> [2D/3D array of float] Install matrix for mode of 'matrix'
            mode -> [str] if 'euler', the classic euler 'ZXZ' rotation transform from RSW to BF is applied
                    if 'ypr', the Yaw(x)–Pitch(z)–Roll(y) rotation transform from RSW to BF is applied
                    if 'matrix', Each column of matrix is the base vector of BF in RSW
                    if 'quaternion', Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format. The quaternion is applied from RSW to BF
            degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees              
        Outputs:
            matrix_trans -> Instance of class MatrixTrans     
        """
        RSW2BF_mat,BF2RSW_mat = frame_trans.RSW_BF_mat(triad,mode,degrees)
        info = {'RSW2BF_mat':RSW2BF_mat,'BF2RSW_mat':BF2RSW_mat,'_transname':'RSW_BF'}
        return MatrixTrans(info)  

    def BF_DF_mat(triad,mode,degrees=True):
        """
        Rotation matrix between the Body-Fixed(BF) reference frame and the Device-Fixed(DF) reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> triad,mode = [20,30,40],'euler' # Classic 'ZXZ' rotation transform from BF to DF
            >>> # triad,mode = [20,30,40],'ypr' # Yaw(x)–Pitch(z)–Roll(y) rotation transform from BF to DF
            >>> # triad,mode = [[0,1,0],[0,0,1],[1,0,0]],'matrix' # Each column of matrix is the base vector of DF in BF
            >>> # triad,mode = [1,2,3,4],'quaternion' # Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format.
            >>> matrix_trans = FrameTrans.BF_DF_mat(triad,mode)
            >>> BF2DF_mat = matrix_trans.BF2DF_mat
            >>> DF2BF_mat = matrix_trans.DF2BF_mat
            >>> print(BF2DF_mat)
            >>> print(DF2BF_mat)
            >>> # For multiple dimension
            >>> from orbdtools import FrameTrans
            >>> triad,mode = [[20,30,40],[60,60,70]],'euler' # Classic 'ZXZ' rotation transform from BF to DF
            >>> #triad,mode = [[20,30,40],[60,60,70]],'ypr' # Yaw(x)–Pitch(z)–Roll(y) rotation transform from BF to DF
            >>> #triad,mode = [[[0,1,0],[0,0,1],[1,0,0]],[[0,-1,0],[0,0,1],[-1,0,0]]],'matrix' # Each column of matrix is the base vector of DF in BF
            >>> #triad,mode = [[1,2,3,4],[5,6,7,8]],'quaternion' # Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format.
            >>> matrix_trans = FrameTrans.BF_DF_mat(triad,mode)
            >>> BF2DF_mat = matrix_trans.BF2DF_mat
            >>> DF2BF_mat = matrix_trans.DF2BF_mat
        Inputs:
            triad -> [1D/2D array of float] Three angles in [rad] or [deg] for mode of 'euler', 'ypr', and 'quaternion'
                  -> [2D/3D array of float] Install matrix for mode of 'matrix'
            mode -> [str] if 'euler', the classic euler 'ZXZ' rotation transform from BF to DF is applied
                    if 'ypr', the Yaw(x)–Pitch(z)–Roll(y) rotation transform from BF to DF is applied
                    if 'matrix', the install matrix of device, Each column of matrix is the base vector of DF in BF
                    if 'quaternion', Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format. The quaternion is applied from BF to DF
            degrees -> [bool,optional,default=True] If True, angles are assumed to be in degrees  
        Outputs:
            matrix_trans -> Instance of class MatrixTrans 
        """
        BF2DF_mat,DF2BF_mat = frame_trans.BF_DF_mat(triad,mode,degrees)
        info = {'BF2DF_mat':BF2DF_mat,'DF2BF_mat':DF2BF_mat,'_transname':'BF_DF'}
        return MatrixTrans(info)

    def ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele,degrees=True):
        """
        Rotation matrix between the ECI reference frame and the Device-Fixed reference frame.

        Usage:
            >>> from orbdtools import FrameTrans
            >>> # For a single payload, a single time
            >>> triad_RSWBF,mode_RSWBF = triad,mode = [20,30,40],'euler' # Classic 'ZXZ' rotation transform from RSW to BF
            >>> triad_BFDF,mode_BFDF = [[0,1,0],[0,0,1],[1,0,0]],'matrix' # Each column of matrix is the base vector of DF in BF 
            >>> orb_ele = [1.5,0.1,20,30,40,50]
            >>> matrix_trans = FrameTrans.ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele)
            >>> ECI2DF_mat = matrix_trans.ECI2DF_mat
            >>> DF2ECI_mat = matrix_trans.DF2ECI_mat
            >>> # For multiple payloads, a single time
            >>> triad_BFDF,mode_BFDF = triad,mode = [[[0,1,0],[0,0,1],[1,0,0]],[[0,-1,0],[0,0,1],[-1,0,0]]],'matrix' # two payload devices
            >>> matrix_trans = FrameTrans.ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele)
            >>> ECI2DF_mat = matrix_trans.ECI2DF_mat
            >>> DF2ECI_mat = matrix_trans.DF2ECI_mat
            >>> # For a single payload, multiple time
            >>> triad_BFDF,mode_BFDF = triad,mode = [[0,1,0],[0,0,1],[1,0,0]],'matrix' 
            >>> orb_ele = [[1.5,0.1,20,30,40,50],[1.6,0.2,20,40,50,60]]
            >>> #triad_RSWBF,mode_RSWBF = triad,mode = [[20,30,40],[60,60,70]],'euler' # Classic 'ZXZ' rotation transform from RSW to BF
            >>> matrix_trans = FrameTrans.ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele)
            >>> ECI2DF_mat = matrix_trans.ECI2DF_mat
            >>> DF2ECI_mat = matrix_trans.DF2ECI_mat
            >>> # For multiple payloads, multiple time
            >>> triad_BFDF,mode_BFDF = triad,mode = [[[0,1,0],[0,0,1],[1,0,0]],[[0,-1,0],[0,0,1],[-1,0,0]]],'matrix' 
            >>> orb_ele = [[1.5,0.1,20,30,40,50],[1.6,0.2,20,40,50,60]]
            >>> #triad_RSWBF,mode_RSWBF = triad,mode = [[20,30,40],[60,60,70]],'euler' # Classic 'ZXZ' rotation transform from RSW to BF
            >>> matrix_trans = FrameTrans.ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele)
            >>> ECI2DF_mat = matrix_trans.ECI2DF_mat
            >>> DF2ECI_mat = matrix_trans.DF2ECI_mat
        Inputs:
            triad_RSWBF/triad_BFDF -> [1D/2D array of float] Three angles in [rad] or [deg] for mode of 'euler', 'ypr' and 'quaternion'
                                   -> [2D/3D array of float] Install matrix for mode of 'matrix' 
            mode_RSWBF/mode_BFDF -> [str] if 'euler', the classic euler 'ZXZ' rotation transform from RSW to BF or from BF to DF is applied
                      if 'ypr', the yaw–pitch–roll rotation transform from RSW to BF or from BF to DF is applied
                      if 'matrix', each column of matrix is the base vector of BF in RSW or DF in BF
                      if 'quaternion', each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format. The quaternion is applied from RSW to BF or from BF to DF
            orb_ele -> [1D/2D array with shape of nx6] Orbit elements in form of [a,ecc,inc,raan,argp,nu], where inc, raan, argp, and nu are in [rad] or [deg]
            degrees -> [bool,optional,default=True] If True, angles are assumed to be in degrees 
        Outputs:
            matrix_trans -> Instance of class MatrixTrans   
        """
        ECI2DF_mat,DF2ECI_mat = frame_trans.ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele,degrees)
        info = {'ECI2DF_mat':ECI2DF_mat,'DF2ECI_mat':DF2ECI_mat,'_transname':'ECI_DF'}
        return MatrixTrans(info)