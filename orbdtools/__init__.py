from .classes.class_arcobs import ArcObs
from .classes.class_tle import TLE
from .utils import data_prepare

# Load and update the EOP file and Leap Second file
data_prepare.iers_load() 