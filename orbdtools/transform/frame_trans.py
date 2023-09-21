import numpy as np
from scipy.spatial.transform import Rotation
from skyfield.framelib import itrs
from skyfield.api import wgs84
from skyfield.sgp4lib import TEME

from ..utils import data_prepare

"""
A brief description of the topocentric horizon coordinate system:
Different reference materials have slightly different definitions to topocentric horizon coordinate system, and special attention should be paid when using them.
For example, NEU(North East Up) coordinate system is used by skyfield and this program, where X points North, Y points East, Z points Up; 
SEZ(South East Zenith) is used by <Bate, Roger R., et al. Fundamentals of astrodynamics. Courier Dover Publications, 2020.> and <Vallado, D. A., and W. D. McClain. "Fundamentals of astrodynamics and applications 4th Edition." (2013).>
ENZ(East North Zenith) is used by <Curtis, Howard. Orbital Mechanics for Engineering Students: Revised Reprint. Butterworth-Heinemann, 2020.>
"""

def reflect_axis(axis):
    """
    Reflect a reference frame by an axis. 
    For example, X Y Z ---> -X Y Z by 'X' axis, X Y Z ---> X -Y Z by 'Y' axis, and X Y Z ---> X Y -Z by 'Z' axis

    Usage:
        >>> M = reflect_axis('X')
        >>> v1 = [1,2,3] # vector in original reference frame
        >>> v2 = M @ v1 # vector in reflected reference frame
        >>> print(v2)
        # we get [-1.,  2.,  3.]
    Inputs:
        axis -> [str] Reflection axis. It can be 'X', 'Y', and 'Z'
    Outputs:
        M -> [2D array] Transformation matrix
    """
    M = np.eye(3)
    if axis == 'X':
        M[0,0] = -1
    elif axis == 'Y':
        M[1,1] = -1
    elif axis == 'Z':
        M[2,2] = -1   
    else:
        raise Exception("reflection axis must be in ['X','Y','Z']")
    return M   

def reflect(seq):
    """
    Reflect a reference frame by a sequence of axes. 

    Usage:
        >>> M = reflect('XXY')
        >>> v1 = [1,2,3] # vector in original reference frame
        >>> v2 = M @ v1 # vector in reflected reference frame
        >>> print(v2)
        # we get [1.,  -2.,  3.]
    Inputs:
        seq -> [str] Sequence of reflection axis, such as 'XXY'
    Outputs:
        M -> [2D array] Transform matrix
    """
    B =  np.eye(3)
    for axis in seq:
        A = reflect_axis(axis)
        B = np.matmul(A,B)
    return B       

def Rot(seq,angles,degrees=True):
    """
    Rotate a reference frame around a sequence of axes.
    Note: The Rot is noly suitable for a right-handed reference frame !!

    Usage:
        >>> seq = 'XY'
        >>> angles = [60,30]
        >>> rotation_matrix = Rot(seq,angles)
    Inputs:
        seq -> [str] Sequence of axes for rotation, such as 'Z' or 'XY'.
        angles -> [float,list of float] Rotation angles in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the rotation angles are assumed to be in degrees
    Outputs:
        rotation_matrix -> [2D array of 3x3] Rotation matrix from the source reference frame to target reference frame   
    """
    if np.isscalar(angles):  
        rotation_matrix = Rotation.from_euler(seq,angles,degrees).as_matrix().T    
    else:  
        angles = np.array(angles)
        if len(seq) > 1:  
            if angles.ndim == 1:
                rotation_matrix = Rotation.from_euler(seq,angles,degrees).as_matrix().T
            elif angles.ndim == 2:
                rotation_matrix = Rotation.from_euler(seq,angles,degrees).as_matrix().transpose(0,2,1)
        else:
            rotation_matrix = Rotation.from_euler(seq,angles,degrees).as_matrix().transpose(0,2,1)
                    
    return rotation_matrix

def euler2vectors(angles,degrees=True):
    """
    Calculate the basis vectors of the source reference frame(XYZ) in target reference frame(xyz). 

    Usage:
        >>> angles = [60,30,40]
        >>> angles = [[60,30,40],[10,120,-70]]
        >>> vector_x,vector_y,vector_z = euler2vectors(angles)
    Inputs:
        angles -> [list of float or 2D array with shape of nx3] Euler angles in sequence of 'ZXZ' in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the euler angles are assumed to be in degrees 
    Outputs:
        vector_x -> [array-like] Basis vector for axis-X
        vector_y -> [array-like] Basis vector for axis-Y
        vector_z -> [array-like] Basis vector for axis-Z
    """
    angles = np.array(angles)
    if len(angles.T) != 3: raise Exception('Euler angles should contain three angles.')
    rotation_matrix = Rotation.from_euler('ZXZ',angles,degrees).as_matrix()

    if angles.ndim == 1:
        vectorX,vectorY,vectorZ = rotation_matrix.T
    elif angles.ndim == 2:
        vectorX,vectorY,vectorZ = rotation_matrix.transpose(2,0,1)

    return vectorX,vectorY,vectorZ    

def lrf_topo_mat(alpha,degrees=True):
    """
    Rotation matrix between the Launch Reference Frame and the Topocentric NEU(North East Up) Reference Frame.

    Usage:
        >>> lrf2topo_mat,topo2lrf_mat = lrf_topo_mat(30)
    Inputs:
        alpha -> [float,list of float] Azimuth of the launch direction in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the azimuth is assumed to be in degrees
    Outputs:
        lrf2topo_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from lrf to topo
        topo2lrf_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from topo to lrf
    """
    if not degrees: 
        alpha = np.rad2deg(alpha)

    # Convert NEZ to right-handed NWZ, then apply rotation
    temp_array1 = np.ones_like(alpha) * 90
    temp_array2 = np.stack([temp_array1,alpha]).T
    topo2lrf_mat =  Rot('XY',temp_array2) @ reflect('Y')

    if temp_array2.ndim == 1:
        lrf2topo_mat = topo2lrf_mat.T
    elif temp_array2.ndim == 2:
        lrf2topo_mat = topo2lrf_mat.transpose(0,2,1)

    return lrf2topo_mat,topo2lrf_mat

def topo_itrf_mat(lon,lat,degrees=True):
    """
    Rotation matrix between the Topocentric NEU(North East Up) and the ITRF(International Terrestrial Reference Frame).

    Usage:
        >>> topo2itrf_mat,itrf2topo_mat = topo_itrf_mat(102.5,25.2)
    Inputs:
        lon -> [float,list of float] Longitude of site in [rad] or [deg]
        lat -> [float,list of float] Latitude of site in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the longitude and latitude are assumed to be in degrees
    Outputs: 
        topo2itrf_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from topo to itrf
        itrf2topo_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from itrf to topo
    """
    if degrees: 
        lon,lat = np.array(lon),np.array(lat) 
    else:
        lon,lat = np.rad2deg(lon),np.rad2deg(lat) 

    loncolat_array = np.stack([lon,90-lat]).T
    itrf2topo_mat = reflect('X') @ Rot('ZY',loncolat_array)

    if loncolat_array.ndim == 1:
        topo2itrf_mat = itrf2topo_mat.T
    elif loncolat_array.ndim == 2:
        topo2itrf_mat = itrf2topo_mat.transpose(0,2,1)

    return topo2itrf_mat,itrf2topo_mat

def gcrf_itrf_mat(ta):
    """
    Rotation matrix between the GCRF(Geocentric Celestial Reference Frame) and the ITRF(International Terrestrial Reference Frame).

    Usage:
        >>> gcrf2itrf_mat,itrf2gcrf_mat = gcrf_itrf_mat(ta)
    Inputs:
        ta -> [array-like of Astropy Time object] Time to make transformation
    Outputs:
        gcrf2itrf_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from gcrf to itrf
        itrf2gcrf_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from itrf to gcrf
    """
    ts = data_prepare.ts
    tf = ts.from_astropy(ta)
        
    if ta.shape:   
        gcrf2itrf_mat = itrs.rotation_at(tf).transpose(2,0,1)
        itrf2gcrf_mat = gcrf2itrf_mat.transpose(0,2,1)
    else:
        gcrf2itrf_mat = itrs.rotation_at(tf)
        itrf2gcrf_mat = gcrf2itrf_mat.T
            
    return gcrf2itrf_mat,itrf2gcrf_mat

def gcrf_topo_mat(lon,lat,ta,degrees=True):
    """
    Rotation matrix between the GCRF(Geocentric Celestial Reference Frame) and the Topocentric NEU(North East Up) Reference Frame.

    Usage:
        >>> gcrf2topo_mat,topo2gcrf_mat = gcrf_topo_mat(lon,lat,ta)
    Inputs:
        lon -> [float,list of float, m=len(lon)] Longitude of site in [rad] or [deg]
        lat -> [float,list of float, m=len(lat)] Latitude of site in [rad] or [deg] 
        ta -> [array-like of Astropy Time object, n=len(ta)] Time to make transformation
        degrees -> [bool,optional,default=True] If True, the longitude and latitude are assumed to be in degrees 
    Outputs:
        gcrf2topo_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from gcrf to topo
        topo2gcrf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from topo to gcrf
    """
    ts = data_prepare.ts
    tf = ts.from_astropy(ta)

    if not degrees:
        lon,lat = np.rad2deg(lon),np.rad2deg(lat) 

    site = wgs84.latlon(lat, lon)

    lonlat_array = np.stack([lon,lat]).T

    if lonlat_array.ndim > 1:
        if ta.shape:
            gcrf2topo_mat,topo2gcrf_mat = [],[]

            for (lon_i,lat_i) in lonlat_array:
                site_i = wgs84.latlon(lat_i, lon_i)
                gcrf2topo_mat_i = site_i.rotation_at(tf).transpose(2,0,1)
                topo2gcrf_mat_i = gcrf2topo_mat_i.transpose(0,2,1)
                gcrf2topo_mat.append(gcrf2topo_mat_i)
                topo2gcrf_mat.append(topo2gcrf_mat_i)
            gcrf2topo_mat,topo2gcrf_mat = np.array(gcrf2topo_mat),np.array(topo2gcrf_mat)    

        else:
            gcrf2topo_mat = site.rotation_at(tf).transpose(2,0,1)
            topo2gcrf_mat = gcrf2topo_mat.transpose(0,2,1)
    else:    
        if ta.shape:
            gcrf2topo_mat = site.rotation_at(tf).transpose(2,0,1)
            topo2gcrf_mat = gcrf2topo_mat.transpose(0,2,1)
        else:     
            gcrf2topo_mat = site.rotation_at(tf)
            topo2gcrf_mat = gcrf2topo_mat.T

    return gcrf2topo_mat,topo2gcrf_mat  

def gcrf_teme_mat(ta):
    """
    Rotation matrix between the GCRF(Geocentric Celestial Reference Frame) and the TEME(True Equator, Mean Equinox) Reference Frame.

    Usage:
        >>> gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta)
    Inputs:
        ta -> [array-like of Astropy Time object] Time to make transformation
    Outputs:
        gcrf2teme_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from gcrf to teme
        teme2gcrf_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from teme to gcrf
    """
    ts = data_prepare.ts
    tf = ts.from_astropy(ta)

    if ta.shape:   
        gcrf2teme_mat = TEME.rotation_at(tf).transpose(2,0,1)
        teme2gcrf_mat = gcrf2teme_mat.transpose(0,2,1)
    else:    
        gcrf2teme_mat = TEME.rotation_at(tf)
        teme2gcrf_mat = gcrf2teme_mat.T

    return gcrf2teme_mat,teme2gcrf_mat

def meme_topo_mat(lon,lat,ta,degrees=True):   
    """
    Rotation matrix between the MEME(Mean Equator, Mean Equinox) Reference Frame and the Topocentric NEU(North East Up) Reference Frame.

    Usage:
        >>> meme2topo_mat,topo2meme_mat = meme_topo_mat(lon,lat,ta)
    Inputs:
        lon -> [float,list of float, m=len(lon)] Longitude of site in [rad] or [deg]
        lat -> [float,list of float, m=len(lat)] Latitude of site in [rad] or [deg]
        ta -> [array-like of Astropy Time object, n=len(ta)] Time to make transformation
        degrees -> [bool,optional,default=True] If True, the longitude and latitude are assumed to be in degrees
    Outputs:
        meme2topo_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from meme to topo
        topo2meme_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from topo to meme
    """
    ts = data_prepare.ts
    tf = ts.from_astropy(ta)

    if not degrees: 
        lon,lat = np.rad2deg(lon),np.rad2deg(lat) 

    gcrf2topo_mat,topo2gcrf_mat = gcrf_topo_mat(lon,lat,ta)
    if ta.shape:
        precession_matrix_t = tf.precession_matrix().transpose(2,0,1)
        if topo2gcrf_mat.ndim > 3:
            topo2meme_mat = precession_matrix_t @ topo2gcrf_mat
            meme2topo_mat = topo2meme_mat.transpose(0,1,3,2)
        else:    
            topo2meme_mat = precession_matrix_t @ topo2gcrf_mat
            meme2topo_mat = topo2meme_mat.transpose(0,2,1)
    else:    
        topo2meme_mat = tf.precession_matrix() @ topo2gcrf_mat

        if topo2meme_mat.ndim > 2 :
            meme2topo_mat = topo2meme_mat.transpose(0,2,1)
        else:
            meme2topo_mat = topo2meme_mat.T
            
    return meme2topo_mat,topo2meme_mat

def lrf_gcrf_mat(lon,lat,alpha,ta,degrees=True):
    """
    Rotation matrix between the Launch Reference Frame and GCRF(Geocentric Celestial Reference Frame).

    Usage:
        >>> lrf2gcrf_mat,gcrf2lrf_mat = lrf_gcrf_mat(lon,lat,alpha,ta)

    Inputs:
        lon -> [float,list of float, m=len(lon)] Longitude of site in [rad] or [deg]
        lat -> [float,list of float, m=len(lat)] Latitude of site in [rad] or [deg]
        alpha -> [float,list of float, m=len(alpha)] Azimuth of the launch direction in [rad] or [deg]
        ta -> [array-like of Astropy Time object, n=len(ta)] Time to make transformation
        degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
    Outputs:
        lrf2gcrf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from lrf to gcrf
        gcrf2lrf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from gcrf to lrf
    """
    if not degrees: 
        lon,lat,alpha = np.rad2deg(lon),np.rad2deg(lat),np.rad2deg(alpha) 

    lonlat_array = np.stack([lon,lat]).T

    lrf2topo_mat,topo2lrf_mat = lrf_topo_mat(alpha)
    gcrf2topo_mat,topo2gcrf_mat = gcrf_topo_mat(lon,lat,ta) 

    if lonlat_array.ndim > 1:
        if ta.shape:
            lrf2gcrf_mat = (topo2gcrf_mat.transpose(1,0,2,3) @ lrf2topo_mat).transpose(1,0,2,3) 
            gcrf2lrf_mat = lrf2gcrf_mat.transpose(0,1,3,2)
        else:    
            lrf2gcrf_mat = topo2gcrf_mat @ lrf2topo_mat
            gcrf2lrf_mat = lrf2gcrf_mat.transpose(0,2,1)
    else:
        lrf2gcrf_mat = topo2gcrf_mat @ lrf2topo_mat

        if ta.shape:
            gcrf2lrf_mat = lrf2gcrf_mat.transpose(0,2,1)
        else:    
            gcrf2lrf_mat = lrf2gcrf_mat.T

    return lrf2gcrf_mat,gcrf2lrf_mat

def lirf_lrf_mat(lon,lat,alpha,ta0,ta,degrees=True):
    """
    Rotation matrix between the Launch Inertial Reference Frame and the Launch Reference Frame.

    Usage:
        >>> lrf2gcrf_mat,gcrf2lrf_mat = lirf_lrf_mat(lon,lat,alpha,ta0,ta)
    Inputs:
        lon -> [float,list of float, m=len(lon)] Longitude of site in [rad] or [deg]
        lat -> [float,list of float, m=len(lat)] Latitude of site in [rad] or [deg]
        alpha -> [float,list of float, m=len(alpha)] Azimuth of the lrf direction in [rad] or [deg]
        ta0 -> [array-like of Astropy Time object, n=len(ta0)] Epoch at which the inertial reference frame defined
        ta -> [array-like of Astropy Time object, n=len(ta)] Time to make transformation
        degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
    Outputs:
        lirf2lrf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from lirf to lrf
        lrf2lirf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from lrf to lirf
    """
    if not degrees: 
        lon,lat,alpha = np.rad2deg(lon),np.rad2deg(lat),np.rad2deg(alpha) 

    lirf2gcrf_mat,gcrf2lirf_mat = lrf_gcrf_mat(lon,lat,alpha,ta0)
    lrf2gcrf_mat,gcrf2lrf_mat = lrf_gcrf_mat(lon,lat,alpha,ta)

    lirf2lrf_mat = gcrf2lrf_mat @ lirf2gcrf_mat

    ndim = lirf2lrf_mat.ndim
    if ndim == 2:
        lrf2lirf_mat = lirf2lrf_mat.T
    elif ndim == 3: 
        lrf2lirf_mat = lirf2lrf_mat.transpose(0,2,1)
    elif ndim == 4:  
        lrf2lirf_mat = lirf2lrf_mat.transpose(0,1,3,2)

    return lirf2lrf_mat,lrf2lirf_mat  

def lrf_meme_mat(lon,lat,alpha,ta,degrees=True):
    """
    Rotation matrix between the Launch Reference Frame and the MEME(Mean Equator, Mean Equinox) Reference Frame.

    Usage:
        >>> lrf2meme_mat,meme2lrf_mat = lrf_meme_mat(lon,lat,alpha,ta)
    Inputs:
        lon -> [float,list of float, m=len(lon)] Longitude of site in [rad] or [deg]
        lat -> [float,list of float, m=len(lat)] Latitude of site in [rad] or [deg]
        alpha -> [float,list of float, m=len(alpha)] Azimuth of the lrf direction in [rad] or [deg]
        ta -> [array-like of Astropy Time object, n=len(ta)] Time to make transformation
        degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
    Outputs:
        lrf2meme_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from lrf to meme
        meme2lrf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from meme to lrf
    """
    if not degrees: 
        lon,lat,alpha = np.rad2deg(lon),np.rad2deg(lat),np.rad2deg(alpha) 

    meme2topo_mat,topo2meme_mat = meme_topo_mat(lon,lat,ta)
    lrf2topo_mat,topo2lrf_mat = lrf_topo_mat(alpha)

    ndim = topo2meme_mat.ndim
    if ndim == 2:
        lrf2meme_mat = topo2meme_mat @ lrf2topo_mat
        meme2lrf_mat = lrf2meme_mat.T
    elif ndim == 3:
        lrf2meme_mat = topo2meme_mat @ lrf2topo_mat
        meme2lrf_mat = lrf2meme_mat.transpose(0,2,1)
    elif ndim == 4:
        lrf2meme_mat = (topo2meme_mat.transpose(1,0,2,3) @ lrf2topo_mat).transpose(1,0,2,3)
        meme2lrf_mat = lrf2meme_mat.transpose(0,1,3,2)
           
    return lrf2meme_mat,meme2lrf_mat

def lrf_teme_mat(lon,lat,alpha,ta,degrees=True):
    """
    Rotation matrix between the Launch Reference Frame and the TEME(True Equator, Mean Equinox) Reference Frame.

    Usage:
        >>> lrf2teme_mat, teme2lrf_mat = lrf_teme_mat(lon,lat,alpha,ta)
    Inputs:
        lon -> [float,list of float, m=len(lon)] Longitude of site in [rad] or [deg]
        lat -> [float,list of float, m=len(lat)] Latitude of site in [rad] or [deg]
        alpha -> [float,list of float, m=len(alpha)] Azimuth of the lrf direction in [rad] or [deg]
        ta -> [array-like of Astropy Time object, n=len(ta)] Time to make transformation
        degrees -> [bool,optional,default=True] If True, the longitude, latitude, and azimuth are assumed to be in degrees
    Outputs:
        lrf2teme_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from lrf to teme
        teme2lrf_mat -> [3D array with shape of 1x3x3 or 4D array with shape of mxnx3x3] Rotation matrix from teme to lrf
    """
    if not degrees: 
        lon,lat,alpha = np.rad2deg(lon),np.rad2deg(lat),np.rad2deg(alpha) 

    lrf2gcrf_mat,gcrf2lrf_mat = lrf_gcrf_mat(lon,lat,alpha,ta)
    gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(ta)
    lrf2teme_mat = gcrf2teme_mat @ lrf2gcrf_mat

    ndim = lrf2gcrf_mat.ndim
    if ndim == 2: 
        teme2lrf_mat = lrf2teme_mat.T
    elif ndim == 3:
        teme2lrf_mat = lrf2teme_mat.transpose(0,2,1)
    elif ndim == 4:
        teme2lrf_mat = lrf2teme_mat.transpose(0,1,3,2)

    return lrf2teme_mat, teme2lrf_mat

def ECI_PQW_mat(inc,raan,argp,degrees=True):  
    """
    Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the PQW(perifocal) reference frame.

    Usage:
        >>> ECI2PQW_mat,PQW2ECI_mat = ECI_PQW_mat(inc,raan,argp)
    Inputs:
        inc -> [float, list of float] Orbital inclination (i) in [rad] or [deg]
        raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
        argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
    Outputs:
        ECI2PQW_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from ECI to PQW
        PQW2ECI_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from PQW to ECI 
    """
    if not degrees:
        inc,raan,argp = np.rad2deg(inc),np.rad2deg(raan),np.rad2deg(argp) 

    temp_array = np.stack([raan,inc,argp]).T
    ECI2PQW_mat = Rot('ZXZ',temp_array)

    if ECI2PQW_mat.ndim == 2:
        PQW2ECI_mat = ECI2PQW_mat.T
    else:
        PQW2ECI_mat = ECI2PQW_mat.transpose(0,2,1)
            
    return ECI2PQW_mat,PQW2ECI_mat

def ECI_RSW_mat(inc,raan,argp,nu,degrees=True):  
    """
    Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the RSW(x:Radial,y:Along-track,z:Cross-track) reference frame.

    Usage:
        >>> ECI2RSW_mat,RSW2ECI_mat = ECI_RSW_mat(inc,raan,argp,nu)
    Inputs:
        inc -> [float, list of float] Orbital inclination (i)  in [rad] or [deg]
        raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
        argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
        nu -> [float, list of float] True anomaly (ν) in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
    Outputs:
        ECI2RSW_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from ECI to RSW
        RSW2ECI_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from RSW to ECI
    """
    if not degrees:
        inc,raan,argp,nu = np.rad2deg(inc),np.rad2deg(raan),np.rad2deg(argp),np.rad2deg(nu)  

    ECI2PQW_mat,PQW2ECI_mat = ECI_PQW_mat(inc,raan,argp)
    ECI2RSW_mat = Rot('Z',nu) @ ECI2PQW_mat

    if ECI2RSW_mat.ndim == 2:
        RSW2ECI_mat = ECI2RSW_mat.T
    else:
        RSW2ECI_mat = ECI2RSW_mat.transpose(0,2,1)
            
    return ECI2RSW_mat,RSW2ECI_mat  

def ECI_NTW_mat(ecc,inc,raan,argp,nu,degrees=True):  
    """
    Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the NTW(x:Normal,y:Tangent,z:Cross-track) reference frame.

    Usage:
        >>> ECI2NTW_mat,NTW2ECI_mat = ECI_NTW_mat(inc,raan,argp,nu)
    Inputs:
        ecc -> [float, list of float] Orbital eccentricity (e)
        inc -> [float, list of float] Orbital inclination (i)  in [rad] or [deg]
        raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
        argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
        nu -> [float, list of float] True anomaly (ν) in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
    Outputs:
        ECI2NTW_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from ECI to NTW
        NTW2ECI_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from NTW to ECI
    """
    if not degrees: 
        inc,raan,argp,nu = np.rad2deg(inc),np.rad2deg(raan),np.rad2deg(argp),np.rad2deg(nu)  

    # Claculate Flight Path Angle
    nu_rad = np.deg2rad(nu)
    fpa = np.arctan2(ecc*np.sin(nu_rad),1+ecc*np.cos(nu_rad))
    ECI2RSW_mat,RSW2ECI_mat = ECI_RSW_mat(inc,raan,argp,nu)
    ECI2NTW_mat = Rot('Z',-fpa,degrees=False) @ ECI2RSW_mat

    if ECI2NTW_mat.ndim == 2:
        NTW2ECI_mat = ECI2NTW_mat.T
    else: 
        NTW2ECI_mat = ECI2NTW_mat.transpose(0,2,1)  

    return ECI2NTW_mat,NTW2ECI_mat  

def ECI_RADAR_mat(inc,raan,argp,nu,degrees=True):  
    """
    Rotation matrix between the Earth Centred Inertial (ECI) reference frame and the RADAR(x:Along-track,y:Cross-track,z:Radial) reference frame.

    Usage:
        >>> ECI2RADAR_mat,RADAR2ECI_mat = ECI_RADAR_mat(inc,raan,argp,nu)
    Inputs:
        inc -> [float, list of float] Orbital inclination (i)  in [rad] or [deg]
        raan -> [float, list of float] Longitude of the ascending node (Ω) in [rad] or [deg]
        argp -> [float, list of float] Argument of periapsis (ω) in [rad] or [deg]
        nu -> [float, list of float] True anomaly (ν) in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees
    Outputs:
        ECI2RADAR_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from ECI to RADAR
        RADAR2ECI_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from RADAR to ECI
    """
    if not degrees:
        inc,raan,argp,nu = np.rad2deg(inc),np.rad2deg(raan),np.rad2deg(argp),np.rad2deg(nu)  

    ECI2RSW_mat,RSW2ECI_mat = ECI_RSW_mat(inc,raan,argp,nu)
    ECI2RADAR_mat = Rot('Z',90) @ Rot('Y',90) @ ECI2RSW_mat

    if ECI2RADAR_mat.ndim == 2:
        RADAR2ECI_mat = ECI2RADAR_mat.T
    else:
        RADAR2ECI_mat = ECI2RADAR_mat.transpose(0,2,1)
            
    return ECI2RADAR_mat,RADAR2ECI_mat

def RSW_BF_mat(triad,mode,degrees=True):  
    """
    Rotation matrix between the RSW(x:Radial,y:Along-track,z:Cross-track) reference frame and the Body-Fixed reference frame.

    Usage:
        >>> RSW2BF_mat,BF2RSW_mat = RSW_BF_mat(triad,mode)
    Inputs:
        triad -> [1D/2D array of float] Three angles in [rad] or [deg] for mode of 'euler', 'ypr' and 'quaternion'
              -> [2D/3D array of float] Install matrix for mode of 'matrix'
        mode -> [str] if 'euler', the classic euler 'ZXZ' rotation transform from RSW to BF is applied
                      if 'ypr', the Yaw(x)–Pitch(z)–Roll(y) rotation transform from RSW to BF is applied
                      if 'matrix', Each column of matrix is the base vector of BF in RSW
                      if 'quaternion', Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format. The quaternion is applied from RSW to BF
        degrees -> [bool,optional,default=True] If True, the angles are assumed to be in degrees              
    Outputs:
        RSW2BF_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from RSW to BF
        BF2RSW_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from BF to RSW 
    """
    triad = np.array(triad)
        
    if not degrees:
        if mode in ['euler','ypr']:
            triad = np.rad2deg(triad)

    if mode == 'euler':
        # alpha,beta,gamma = triad
        RSW2BF_mat = Rot('ZXZ',triad)
    elif mode == 'ypr':
        RSW2BF_mat = Rot('XZY',triad)
    elif mode == 'matrix':
        triad_ndim = triad.ndim
        if triad_ndim == 2:
            RSW2BF_mat = triad.T
        elif triad_ndim == 3:
            RSW2BF_mat = triad.transpose(0,2,1) 
    elif mode == 'quaternion':  
        RSW2BF_mat = Rotation.from_quat(triad).as_matrix()     
    else:
        raise Exception("'mode' must be in ['euler','ypr','matrix','quatern']")     

    RSW2BF_mat_ndim = RSW2BF_mat.ndim
    if RSW2BF_mat_ndim == 2:    
        BF2RSW_mat = RSW2BF_mat.T   
    elif RSW2BF_mat_ndim == 3:
        BF2RSW_mat = RSW2BF_mat.transpose(0,2,1)  
             
    return RSW2BF_mat,BF2RSW_mat 

def BF_DF_mat(triad,mode,degrees=True):
    """
    Rotation matrix between the Body-Fixed(BF) reference frame and the Device-Fixed(DF) reference frame.

    Usage:
        >>> BF2DF_mat,DF2BF_mat = BF_DF_mat(triad,mode)
    Inputs:
        triad -> [1D/2D array of float] Three angles in [rad] or [deg] for mode of 'euler', 'ypr' and 'quaternion'
              -> [2D/3D array of float] Install matrix for mode of 'matrix' 
        mode -> [str] if 'euler', the classic euler 'ZXZ' rotation transform from BF to DF is applied
                      if 'ypr', the Yaw(x)–Pitch(z)–Roll(y) rotation transform from BF to DF is applied
                      if 'matrix', the install matrix of device, Each column of matrix is the base vector of DF in BF
                      if 'quaternion', Each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format. The quaternion is applied from BF to DF
        degrees -> [bool,optional,default=True] If True, angles are assumed to be in degrees  
    Outputs:
        BF2DF_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from BF to DF
        DF2BF_mat -> [2D array or 3D array with shape of nx3x3] Rotation matrix from DF to BF
    """
    BF2DF_mat,DF2BF_mat = RSW_BF_mat(triad,mode,degrees)
             
    return BF2DF_mat,DF2BF_mat 

def ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele,degrees=True):
    """
    Rotation matrix between the ECI reference frame and the Device-Fixed reference frame.

    Usage:
        >>> ECI2DF_mat,DF2ECI_mat = ECI_DF_mat(triad_RSWBF,mode_RSWBF,triad_BFDF,mode_BFDF,orb_ele)
    Inputs:
        triad_RSWBF/triad_BFDF -> [1D/2D array of float] Three angles in [rad] or [deg] for mode of 'euler', 'ypr' and 'quaternion'
                               -> [2D/3D array of float] Install matrix for mode of 'matrix' 
        mode_RSWBF/mode_BFDF -> [str] if 'euler', the classic euler 'ZXZ' rotation transform from RSW to BF or from BF to DF is applied
                      if 'ypr', the yaw–pitch–roll rotation transform from RSW to BF or from BF to DF is applied
                      if 'matrix', each column of matrix is the base vector of BF in RSW or DF in BF
                      if 'quaternion', each row is a (possibly non-unit norm) quaternion in scalar-last (x, y, z, w) format. The quaternion is applied from RSW to BF or from BF to DF
        orb_ele ->  [1D/2D array with shape of nx6] Orbit elements in form of [a,ecc,inc,raan,argp,nu], where inc, raan, argp, and nu are in [rad] or [deg]
        degrees -> [bool,optional,default=True] If True, angles are assumed to be in degrees 
    Outputs:
        ECI2DF_mat -> [4D array with shape of nx3x3] Rotation matrix from ECI to DF with shape of mxnx3x3
        DF2ECI_mat -> [4D array with shape of nx3x3] Rotation matrix from DF to ECI with shape of mxnx3x3 
    """
    triad_RSWBF = np.array(triad_RSWBF)
    triad_BFDF = np.array(triad_BFDF)
    orb_ele = np.array(orb_ele)
    inc,raan,argp,nu = orb_ele.T[2:]
        
    if not degrees:
        if mode_RSWBF in ['euler','ypr']:
            triad_RSWBF = np.rad2deg(triad_RSWBF)
        if mode_BFDF == 'euler':
            triad_BFDF = np.rad2deg(triad_BFDF)

        inc,raan,argp,nu = np.rad2deg(orb_ele).T[2:]

    ECI2RSW_mat,RSW2ECI_mat = ECI_RSW_mat(inc,raan,argp,nu) # nx3x3
    RSW2BF_mat,BF2RSW_mat = RSW_BF_mat(triad_RSWBF,mode_RSWBF) # nx3x3
    BF2DF_mat, DF2BF_mat = BF_DF_mat(triad_BFDF,mode_BFDF) # mx3x3

    BF2DF_mat_ndim = BF2DF_mat.ndim
    RSW2BF_mat_ndim = RSW2BF_mat.ndim

    if BF2DF_mat_ndim == 2: BF2DF_mat = np.array([BF2DF_mat]) # 1x3x3
    if RSW2BF_mat_ndim == 2:
        ECI2DF_mat = (BF2DF_mat @ RSW2BF_mat)[:,None] @ ECI2RSW_mat
    elif RSW2BF_mat_ndim == 3:
        ECI2DF_mat = (BF2DF_mat[:,None] @ RSW2BF_mat) @ ECI2RSW_mat  

    DF2ECI_mat = ECI2DF_mat.transpose(0,1,3,2) 
        
    return ECI2DF_mat,DF2ECI_mat      