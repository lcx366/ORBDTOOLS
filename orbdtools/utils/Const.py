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

# non-dimensional quantities
L_nd = Re
# M_nd = Me
T_nd = np.sqrt(Re**3/mu)
v_nd = L_nd / T_nd # in [L]/[T]
G_nd = 1 # in [L]^3/([M]*[T]^2)
mu_nd = 1 # in [L]^3/[T]^2

au = 1.495978707e8 # Astronomical unit, [km]
S_R = 4.56e-6 # solar radiation pressure, [N/m^2]
R_sun = 6.957e5 # Average radius of sun in [km]
R_moon = 1737.4 # Average radius of moon in [km]
stellar_year = 365.2564 # mean solar day
n_sunsync = 1.99096871e-7 # the precession rate of Sun-synchronous orbit in rad/s

# Standard gravitational parameter
mu_sun = 1.32712440018e11 # GM for Sun, [km^3/s^2]
mu_moon = 4904.8695 # GM for Moon, [km^3/s^2]
mu_jup = 1.26686534e8 # GM for Jupiter, [km^3/s^2]

# Time 
t_J2000 = Time('2000-1-1 12:00:00',scale='ut1')
jd_J2000 = t_J2000.jd

# other constants
C_D = 2.2 # drag codfficient
C_R = 1.5 # radiation pressure codfficient between 1 and 2; 
          # 1 for blackbody, which means absorbing all of the momentum of the incident photon stream
          # 2 means reflecting all of the momentum of the incident photon stream

# parameters configuration for targets
m = 100 # mass, [kg]
A_D = np.pi/4 # section area for drag, [m^2]
A_R = np.pi/4 # section area for solar radiation, [m^2]
