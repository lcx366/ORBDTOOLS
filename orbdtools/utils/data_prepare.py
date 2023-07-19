from skyfield.api import Loader,load
from skyfield.data import iers as iers_skyfield
from astropy.utils import iers as iers_astropy

from .data_download import download_iers

def iers_load():
    """
    Load the EOP and Leap Second file from IERS for time system transformation and coordinate system transformation.
    """

    global ts

    # load the EOP file
    dir_iers,eop_file,leapsecond_file = download_iers()
    load_iers = Loader(dir_iers)
    ts = load_iers.timescale(builtin=False) 
    with load.open(eop_file) as f:
        pm_data = iers_skyfield.parse_x_y_dut1_from_finals_all(f)
    iers_skyfield.install_polar_motion_table(ts, pm_data)  

    # for astropy
    iers_astropy.conf.auto_download = False
    iers_a = iers_astropy.IERS_A.open(eop_file)
    leapsecond = iers_astropy.LeapSeconds.from_iers_leap_seconds(leapsecond_file)
    eop_table = iers_astropy.earth_orientation_table.set(iers_a)

