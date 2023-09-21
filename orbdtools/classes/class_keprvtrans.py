from ..transform import kep_rv_trans                                                

class KeprvTrans(object):
    """
    Class of transformation between classical orbital elements and state vectors.
    """

    def coe2rv(coe,mu,degrees=True): 
        """
        Transform classical orbital elements to state vectors.

        Usage:
            >>> from orbdtools import KeprvTrans
            >>> import numpy as np
            >>> # classical orbital elements in form of [a, e, i, Ω, ω, ν]
            >>> coe = np.array([7000,0.01,50,100,30,210]) # semi major axis is in [km], and angles are in [deg]
            >>> mu = 398600.4418 # GM for the Reference Earth Model - WGS84, [km^3/s^2] 
            >>> rv = KeprvTrans.coe2rv(coe,mu)
            >>> print(rv)
            >>> # [ 4.48364689e+03 -2.79409408e+03 -4.68399786e+03  1.21885031e+00 6.81282168e+00 -2.84038655e+00]
            >>> # For non-dimensional/dimensionless orbital elements
            >>> coe_nd = np.array([1.0974,0.01,50,100,30,210]) # semi major axis is in non-dimensional length unit [L_nd]. For the Reference Earth Model - WGS84, [L_nd]=6378.137 [km] as the equatorial radius.  
            >>> mu_nd = 1.5 # GM of the central attraction in non-dimensional unit [mu_nd]. For the Reference Earth Model - WGS84, [mu_nd]=398600.4418 [km^3/s^2].
            >>> rv_nd = KeprvTrans.coe2rv(coe_nd,mu_nd)
            >>> print(rv_nd)
            >>> # [ 0.70290773 -0.43803412 -0.73431704  0.18883985  1.05552933 -0.44006895]
        Inputs:
            coe -> [array-like] classical orbital elements in form of [a, e, i, Ω, ω, ν], where
                a: semi major axis
                e: eccentricity
                i: inclination, [rad] or [deg]
                Ω: right ascension of the ascending node, [rad] or [deg]
                ω: argument of perigee, [rad] or [deg]
                v: true anomaly, [rad] or [deg]
            mu -> [float] GM of the central attraction
            degrees -> [bool,optional,default=True] unit of i, Ω, ω, ν
        Outputs:
            rv -> [array-like] State vector in form of [x, y, z, vx, vy, vz] 
        """
        rv = kep_rv_trans.coe2rv(coe,mu,degrees)
        return rv

    def rv2coe(rvs,mu,degrees=True,tol=1e-9):
        """
        Transform state vectors to classical orbital elements.

        Usage:
            >>> from orbdtools import KeprvTrans
            >>> import numpy as np
            >>> rvs = np.array([[ 4.48e+03, -2.79e+03, -4.68e+03,  1.22e+00,6.81e+00, -2.84e+00],[ 5.48e+03, -3.79e+03, -5.68e+03,  1.52e+00,7.81e+00, -3.84e+00]])
            >>> mu = 398600.4418 # GM for the Reference Earth Model - WGS84, [km^3/s^2] 
            >>> coe = KeprvTrans.rv2coe(rvs,mu)
            >>> print(coe)
            >>> # [[6.98242989e+03 1.12195221e-02 4.99946032e+01 9.99954947e+01 3.60288078e+01 2.03986947e+02] [3.06604952e+04 7.14457354e-01 5.11123776e+01 1.01894775e+02 2.35493764e+02 9.61451955e-01]]
            >>> # For non-dimensional/dimensionless state vectors
            >>> # The non-dimensional time unit [T_nd] is defined by sqrt([L_nd]**3/[mu_nd]), and non-dimensional velocity unit [v_nd] is defined by [L_nd]/[T_nd]
            >>> rvs_nd = np.array([[ 0.70239946, -0.43743181, -0.73375658,  0.15432556,  0.86144022,-0.35924967],[ 0.85918506, -0.5942174 , -0.89054218,  0.19227447,  0.98793658,-0.48574603]])
            >>> mu_nd = 1.5 
            >>> coe_nd = KeprvTrans.rv2coe(rvs_nd,mu_nd)
            >>> print(coe_nd)
        Inputs:
            rv -> [array-like,float] state vector
            mu -> [float] GM of the central attraction
            degrees -> [bool,optional,default=True] unit of i, Ω, ω, ν  
            tol -> [float,optional,default=1e-9] Threshold for small eccentricity or small inclination. 
                   If the eccentricity is less than the threshold, it is treated as a circular orbit;
                   if the orbital inclination is less than the threshold, it is treated as a equatorial orbit.      
        Outputs:
            coe -> [array-like,float] classical orbital elements, where
                a: semi major axis
                e: eccentricity
                i: inclination, [deg] or [rad]
                Ω: right ascension of the ascending node, [deg] or [rad]
                ω: argument of perigee, [deg] or [rad]
                v: true anomaly, [deg] or [rad]
        """  
        coe = kep_rv_trans.rv2coe(rvs,mu,degrees,tol) 

        return coe 