import numpy as np
from astropy.time import Time

# basic physical constants for the Reference Earth Model - WGS84
mu = 398600.4418 # GM, [km^3/s^2]  
Re = 6378.137 # Equatorial radius of the Earth from WGS84 ellipsoid model, [km]
f = 1/298.257223563 # flattening
Re_V = 6371.0008 # volumetric radius, [km] 
rot = 7.292115e-5 # Rotational speed of the Earth from WGS84 ellipsoid model, [rad/s]
J2 = 1.08262982e-3

# basic physical constants for the Reference Earth Model - WGS72
mu_sgp4 = 398600.8 # GM, [km^3/s^2] 
Re_sgp4 = 6378.135 # Equatorial radius, [km]
T_sgp4 = np.sqrt(Re_sgp4**3/mu_sgp4) # [s]
v_sgp4 = Re_sgp4/T_sgp4

ke_sgp4 = 1
k2_sgp4 = 5.41308e-4

# non-dimensional units
L_nd_unit = Re
mu_nd_unit = mu # in [L_nd]^3/[T_nd]^2
T_nd_unit = np.sqrt(L_nd_unit**3/mu_nd_unit) # in ([L_nd]^3/[mu_nd])**0.5
v_nd_unit = L_nd_unit / T_nd_unit # in [L_nd]/[T_nd]
G_nd_unit = 1 # in [L_nd]^3/([M_nd]*[T_nd]^2)

au = 1.495978707e8 # Astronomical unit, [km]
S_R = 4.56e-6 # solar radiation pressure, [N/m^2]
stellar_year = 365.2564 # mean solar day
n_sunsync = 1.99096871e-7 # the precession rate of Sun-synchronous orbit in rad/s


# Time 
t_J2000 = Time('2000-1-1 12:00:00',scale='ut1')
jd_J2000 = t_J2000.jd

# other constants
twopi = 2 * np.pi
