from sgp4.api import Satrec
from skyfield.api import EarthSatellite
import numpy as np

from ..utils import data_prepare

def load_tle_file(filename):
    """
    Load and parse a TLE file.

    Usage:
        >>> sats_Satrec,sats_EarthSatellite = load_tle_file('test/tle_20220524.txt')
    Inputs:
        filename -> [str] path of the TLE file
    Outputs:
        sats_Satrec -> array of objects of class Satrec
        sats_EarthSatellite -> array of objects of class EarthSatellite
    """
    with open(filename,'rb') as f:
        sats_list = list(parse_tle_file(f))  
    sats_Satrec,sats_EarthSatellite = np.array(sats_list).T    
    return sats_Satrec,sats_EarthSatellite
    
def parse_tle_file(lines):
    """
    Parse lines of TLE data, yielding a sequence of satellites.
    """
    ts = data_prepare.ts
    b0 = b1 = b''
    for b2 in lines:
        if (b2.startswith(b'2 ') and len(b2) >= 69 and b1.startswith(b'1 ') and len(b1) >= 69):
            if b0:
                b0 = b0.rstrip(b' \n\r')
                if b0.startswith(b'0 '): b0 = b0[2:]  # Spacetrack 3-line format
            line1 = b1.decode('ascii')
            line2 = b2.decode('ascii')
            sat_Satrec = Satrec.twoline2rv(line1, line2)
            sats_EarthSatellite = EarthSatellite.from_satrec(sat_Satrec, ts)
            yield sat_Satrec,sats_EarthSatellite
            b0 = b1 = b''
        else:
            b0,b1 = b1,b2
    