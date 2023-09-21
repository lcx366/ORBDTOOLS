import numpy as np
from scipy.optimize import fsolve
from sgp4.api import Satrec,WGS72
from sgp4 import exporter
 
from ...utils import Const

def a0_del0(a0,a0dp,ecc,inc):
    """
    from sgp4 
    """
    k2 = Const.k2_sgp4
    
    cosi0 = np.cos(inc)
    x3thm1 = 3 * cosi0**2 - 1
    del0 = 1.5 * k2 * x3thm1 / a0**2 / (1 - ecc**2)**1.5
    return a0 / (1 - del0) - a0dp

def a1_del1(a1,a0,ecc,inc):
    """
    from sgp4 
    """
    k2 = Const.k2_sgp4
    
    cosi0 = np.cos(inc)
    x3thm1 = 3 * cosi0**2 - 1
    del1 = 1.5 * k2 * x3thm1 / a1**2 / (1 - ecc**2)**1.5
    
    return a1*(1 - del1/3 - del1**2 - 134/81*del1**3) - a0

def no_kozai_calculate(am,eccm,incm):
    """
    calculate no_kozai in sgp4 initialization
    am : secular a
    eccm: secular ecc
    incm: secular inc
    """
    
    ke = Const.ke_sgp4
    
    a0 = fsolve(a0_del0,am,args=(am,eccm,incm))
    a1 = fsolve(a1_del1,a0,args=(a0,eccm,incm))
    
    no_kozai = np.sqrt(1/a1**3)/Const.T_sgp4*60 # in [rad/min]

    if len(no_kozai) > 1:
        raise Exception("only one solution is valid")
    else:
        no_kozai = no_kozai[0]    
    
    return no_kozai 

def sgp4init(a,ecc,inc,raan,argp,M,epoch,bstar=0,ndot=0,nddot=0,satnum=-1):
    
    """
    a: major axis in nd
    epoch: UTC
    """
    
    epochdays = epoch.mjd - 33281.0 # 33281.0 is the mjd of 1949 December 31 00:00 UT1

    no_kozai = no_kozai_calculate(a,ecc,inc)
    
    satrec = Satrec()
    satrec.sgp4init(
    WGS72,           # gravity model
    'i',             # 'a' = old AFSPC mode, 'i' = improved mode
    satnum,          # satnum: Satellite number
    epochdays,       # epoch: days since 1949 December 31 00:00 UT
    bstar,           # bstar: drag coefficient (/earth radii)
    ndot,            # ndot: ballistic coefficient or one half the first time derivative of the mean motion (radians/minute^2)
    nddot,           # nddot: one sixth the second derivative of mean motion (radians/minute^3)
    ecc,             # ecco: eccentricity
    argp,            # argpo: argument of perigee (radians)
    inc,             # inclo: inclination (radians)
    M,               # mo: mean anomaly (radians)
    no_kozai,        # no_kozai: mean motion (radians/minute)
    raan,            # nodeo: right ascension of ascending node (radians)
    )
    
    return satrec    


def sgp4init_no_kozai(no_kozai,ecc,inc,raan,argp,M,epoch,bstar=0,ndot=0,nddot=0,satnum=-1):
    
    """
    a: major axis in nd
    epoch: UTC
    """
    
    epochdays = epoch.mjd - 33281.0 # 33281.0 is the mjd of 1949 December 31 00:00 UT1
    
    satrec = Satrec()
    satrec.sgp4init(
    WGS72,           # gravity model
    'i',             # 'a' = old AFSPC mode, 'i' = improved mode
    satnum,          # satnum: Satellite number
    epochdays,       # epoch: days since 1949 December 31 00:00 UT
    bstar,           # bstar: drag coefficient (/earth radii)
    ndot,            # ndot: ballistic coefficient or one half the first time derivative of the mean motion (radians/minute^2)
    nddot,           # nddot: one sixth the second derivative of mean motion (radians/minute^3)
    ecc,             # ecco: eccentricity
    argp,            # argpo: argument of perigee (radians)
    inc,             # inclo: inclination (radians)
    M,               # mo: mean anomaly (radians)
    no_kozai,        # no_kozai: mean motion (radians/minute)
    raan,            # nodeo: right ascension of ascending node (radians)
    )
    
    return satrec     

def export_tle(no_kozai,ecc,inc,raan,argp,M,epoch,bstar=0,ndot=0,nddot=0,satnum=-1,classification='U',intldesg='YYXXXA',elnum=1,revnum=0):
    satrec = sgp4init_no_kozai(no_kozai,ecc,inc,raan,argp,M,epoch,bstar,ndot,nddot,satnum)
    satrec.classification = classification
    satrec.elnum = elnum
    satrec.revnum = revnum
    satrec.intldesg = intldesg

    line1, line2 = exporter.export_tle(satrec)   

    return [line1,line2]