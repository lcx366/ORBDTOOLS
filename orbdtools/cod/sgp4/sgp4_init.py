import numpy as np
from scipy.optimize import fsolve
from sgp4.api import Satrec,WGS72
from sgp4 import exporter
 
from ...utils import Const

def a0_del0(a0,a0dp,ecc,inc):
    """
    Solve for a0 from mean semimajor a0'' in Brouwer theory.

    Reference:
        Hoots F R. Spacetrack report no. 3, models for propagation of norad element sets[J]. http://www. itc. nl/-bakker/orbit. html, 1980.    
    """
    k2 = Const.k2_sgp4
    
    cosi0 = np.cos(inc)
    x3thm1 = 3 * cosi0**2 - 1
    del0 = 1.5 * k2 * x3thm1 / a0**2 / (1 - ecc**2)**1.5
    return a0 / (1 - del0) - a0dp

def a1_del1(a1,a0,ecc,inc):
    """
    Solve for a1 from a0.

    Reference:
        Hoots F R. Spacetrack report no. 3, models for propagation of norad element sets[J]. http://www. itc. nl/-bakker/orbit. html, 1980.    
    """
    k2 = Const.k2_sgp4
    
    cosi0 = np.cos(inc)
    x3thm1 = 3 * cosi0**2 - 1
    del1 = 1.5 * k2 * x3thm1 / a1**2 / (1 - ecc**2)**1.5
    
    return a1*(1 - del1/3 - del1**2 - 134/81*del1**3) - a0

def no_kozai_calculate(am,eccm,incm):
    """
    Calculate initial mean motion in Kozai theory.

    Usage:
        >>> no_kozai = no_kozai_calculate(am,eccm,incm)
    Inputs:
        am -> [float] Initial mean semimajor, equivalent to a0dp
        eccm -> [float] Initial mean eccentricity, equivalent to e0
        incm -> [float] Initial mean inclination, equivalent to i0
    Outputs:
        no_kozai -> [float] Initial mean motion
    Reference:
        Hoots F R. Spacetrack report no. 3, models for propagation of norad element sets[J]. http://www. itc. nl/-bakker/orbit. html, 1980.    
    """
    ke = Const.ke_sgp4
    
    a0 = fsolve(a0_del0,am,args=(am,eccm,incm))
    a1 = fsolve(a1_del1,a0,args=(a0,eccm,incm))
    
    no_kozai = np.sqrt(1/a1**3)/Const.T_sgp4*60 # in [rad/min]

    if len(no_kozai) > 1:
        raise Exception("Multiple solutions exist, but only one is valid.")
    else:
        no_kozai = no_kozai[0]    
    
    return no_kozai 

def sgp4init(a,ecc,inc,raan,argp,M,epoch,bstar=0,ndot=0,nddot=0,satnum=-1):
    """
    Initialize SGP4 using specified orbital elements in TEME.

    Usage:
        >>> satrec = sgp4init(a,ecc,inc,raan,argp,M,epoch)
    Inputs:
        a -> [float] Initial mean Semimajor in [Length unit], 
        where [Length unit] is defined as the equatorial radius of the earth from WGS72, i.e., 6378.135km. 
        ecc -> [float] Initial mean Eccentricity
        inc -> [float] Initial mean Inclination in [rad]
        raan -> [float] Initial mean Right Ascension of Ascending Node in [rad]
        argp -> [float] Initial mean Argument of Perigee in [rad]
        M -> [float] Initial mean Mean Anomaly in [rad]
        epoch -> [Astropy Time] Epoch of the orbital elements
        bstar -> [float,optional,default=0] The drag term, or radiation pressure coefficient, in unit of [1/earth radii]
        ndot -> [float,optional,default=0] Ballistic coefficient or one half the first time derivative of the mean motion in unit of [rad/min^2]
        nddot -> [float,optional,default=0] One sixth the second derivative of mean motion in unit of [rad/min^3]
        satnum -> [int,optional,default=-1] Satellite number
    Outputs:
        satrec -> object of class Satrec
    """
    epochdays = epoch.mjd - 33281.0 # Days since 1949 December 31 00:00 UT

    # Calculate the initial mean motion in Kozai theory.
    no_kozai = no_kozai_calculate(a,ecc,inc) # in [rad/min]
    
    satrec = Satrec()
    satrec.sgp4init(
    WGS72,           # gravity model
    'i',             # 'a' = old AFSPC mode, 'i' = improved mode
    satnum,epochdays,bstar,ndot,nddot,ecc,argp,inc,M,no_kozai,raan)
    
    return satrec    

def sgp4init_no_kozai(no_kozai,ecc,inc,raan,argp,M,epoch,bstar=0,ndot=0,nddot=0,satnum=-1):
    """
    Initialize SGP4 using specified orbital elements in TEME.

    Usage:
        >>> satrec = sgp4init(a,ecc,inc,raan,argp,M,epoch)
    Inputs:
        no_kozai -> [float] Initial mean motion in unit of [rad/min]
        ecc -> [float] Initial mean Eccentricity
        inc -> [float] Initial mean Inclination in [rad]
        raan -> [float] Initial mean Right Ascension of Ascending Node in [rad]
        argp -> [float] Initial mean Argument of Perigee in [rad]
        M -> [float] Initial mean Mean Anomaly in [rad]
        epoch -> [Astropy Time] Epoch of the orbital elements
        bstar -> [float,optional,default=0] The SGP4 type drag coefficient, in unit of [1/earth radii]
        ndot -> [float,optional,default=0] Ballistic coefficient or one half the first time derivative of the mean motion in unit of [rad/min^2]
        nddot -> [float,optional,default=0] One sixth the second derivative of mean motion in unit of [rad/min^3]
        satnum -> [int,optional,default=-1] Satellite number
    Outputs:
        satrec -> object of class Satrec
    """
    epochdays = epoch.mjd - 33281.0 # Days since 1949 December 31 00:00 UT
    
    satrec = Satrec()
    satrec.sgp4init(
    WGS72,           # gravity model
    'i',             # 'a' = old AFSPC mode, 'i' = improved mode
    satnum,epochdays,bstar,ndot,nddot,ecc,argp,inc,M,no_kozai,raan)
    
    return satrec     

def export_tle(no_kozai,ecc,inc,raan,argp,M,epoch,bstar=0,ndot=0,nddot=0,satnum=-1,classification='U',intldesg='YYXXXA',elnum=1,revnum=1000):
    """
    Export orbital elements in form of TLE.

    Usage:
        tle_str = export_tle(no_kozai,ecc,inc,raan,argp,M,epoch)
    Inputs:
        no_kozai -> [float] Initial mean motion in unit of [rad/min]
        ecc -> [float] Initial mean Eccentricity
        inc -> [float] Initial mean Inclination in [rad]
        raan -> [float] Initial mean Right Ascension of Ascending Node in [rad]
        argp -> [float] Initial mean Argument of Perigee in [rad]
        M -> [float] Initial mean Mean Anomaly in [rad]
        epoch -> [Astropy Time] Epoch of the orbital elements
        bstar -> [float,optional,default=0] The SGP4 type drag coefficient, in unit of [1/earth radii]
        ndot -> [float,optional,default=0] Ballistic coefficient or one half the first time derivative of the mean motion in unit of [rad/min^2]
        nddot -> [float,optional,default=0] One sixth the second derivative of mean motion in unit of [rad/min^3]
        satnum -> [int,optional,default=-1] Satellite number
        classification -> [str,optional,default='U'] Classification of the space object. Avaliable options are 'U', 'C', and 'S', which indicate Unclassified, Classified, and Secret respectively.
        intldesg -> [str,optional,default='YYXXXA'] International Designator, where YY is the last two digits of launch year, XXX is the launch number of the year, and A is the piece of the launch.
        elnum -> [int,optional,default=0] Element set number. Incremented when a new TLE is generated for this object.
        revnum -> [int,optional,default=1000] Revolution number at epoch
    Outputs:
        tle_str -> [tuple of str] TLE    
    """
    satrec = sgp4init_no_kozai(no_kozai,ecc,inc,raan,argp,M,epoch,bstar,ndot,nddot,satnum)
    satrec.classification = classification
    satrec.elnum = elnum
    satrec.revnum = revnum
    satrec.intldesg = intldesg

    tle_str = exporter.export_tle(satrec)   

    return tle_str