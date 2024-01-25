from ..transform import kep_rv_trans                                              

class KeprvTrans(object):
    """
    Class of transformation between orbital elements and state vectors.
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
            >>> # =======
            >>> # For non-dimensional/dimensionless orbital elements
            >>> coe_nd = np.array([1.0974,0.01,50,100,30,210]) # semi major axis is in non-dimensional length unit [L_nd]. For the Reference Earth Model - WGS84, [L_nd]=6378.137 [km] as the equatorial radius.  
            >>> mu_nd = 1.5 # GM of the central attraction in non-dimensional unit [mu_nd]. For the Reference Earth Model - WGS84, [mu_nd]=398600.4418 [km^3/s^2].
            >>> rv_nd = KeprvTrans.coe2rv(coe_nd,mu_nd)
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
            >>> # =======
            >>> # For non-dimensional/dimensionless state vectors
            >>> # The non-dimensional time unit [T_nd] is defined by sqrt([L_nd]**3/[mu_nd]), and non-dimensional velocity unit [v_nd] is defined by [L_nd]/[T_nd]
            >>> rvs_nd = np.array([[ 0.70239946, -0.43743181, -0.73375658,  0.15432556,  0.86144022,-0.35924967],[ 0.85918506, -0.5942174 , -0.89054218,  0.19227447,  0.98793658,-0.48574603]])
            >>> mu_nd = 1.5 
            >>> coe_nd = KeprvTrans.rv2coe(rvs_nd,mu_nd)
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

    def mee2rv(mee,mu,degrees=True):
        """
        Transform modified equinoctial elements to state vectors.

        Usage:
            >>> from orbdtools import KeprvTrans
            >>> import numpy as np
            >>> # modified equinoctial elements in form of [p, f, g, h, k, L]
            >>> mee = np.array([6999.3,-6.43e-3,7.66e-3,8.10e-2,0.46,340]) # p and L are in [km] and [degrees] respectively, and the units of other parameters are dimensionless.
            >>> mu = 398600.4418 # GM for the Reference Earth Model - WGS84, [km^3/s^2] 
            >>> rv = KeprvTrans.mee2rv(mee,mu)
            >>> # =======
            >>> # For non-dimensional/dimensionless orbital elements
            >>> mee_nd = np.array([1.097,-6.43e-3,7.66e-3,8.10e-2,0.46,340]) # Semi-latus rectum is in non-dimensional length unit [L_nd]. For the Reference Earth Model - WGS84, [L_nd]=6378.137 [km] as the equatorial radius.  
            >>> mu_nd = 1.5 # GM of the central attraction in non-dimensional unit [mu_nd]. For the Reference Earth Model - WGS84, [mu_nd]=398600.4418 [km^3/s^2].
            >>> rv_nd = KeprvTrans.mee2rv(mee_nd,mu_nd)
        Inputs:
            mee -> [array-like] modified equinoctial elements in form of [p, f, g, h, k, L], where
                p -> [array-like,float] Semi-latus rectum, p = a*(1-e**2)
                f -> [array-like,float] x components of the eccentricity vector in the orbital frame, f = e * cos(Ω + ω)
                g -> [array-like,float] y components of the eccentricity vector in the orbital frame, g = e * sin(Ω + ω)
                h -> [array-like,float] x components of the node vector in the orbital frame, h = tan(i/2) * cos(Ω)
                k -> [array-like,float] y components of the node vector in the orbital frame, k = tan(i/2) * sin(Ω)
                L -> [array-like,float] True Longitude, [radians] or [deg], L = Ω + ω + ν
            mu -> [float] GM of the central attraction
            degrees -> [bool,optional,default=True] units of L
        Outputs:
            rv -> [array-like] state vector in form of [x, y, z, vx, vy, vz]
        """    
        rv = kep_rv_trans.mee2rv(mee,mu,degrees)
        return rv

    def rv2mee(rvs,mu,degrees=True):
        """
        Transform state vectors to modified equinoctial elements.

        Usage:
            >>> from orbdtools import KeprvTrans
            >>> import numpy as np
            >>> rvs = np.array([[ 4.48e+03, -2.79e+03, -4.68e+03,  1.22e+00,6.81e+00, -2.84e+00],[ 5.48e+03, -3.79e+03, -5.68e+03,  1.52e+00,7.81e+00, -3.84e+00]])
            >>> mu = 398600.4418 # GM for the Reference Earth Model - WGS84, [km^3/s^2] 
            >>> mee = KeprvTrans.rv2mee(rvs,mu)
            >>> # =======
            >>> # For non-dimensional/dimensionless state vectors
            >>> # The non-dimensional time unit [T_nd] is defined by sqrt([L_nd]**3/[mu_nd]), and non-dimensional velocity unit [v_nd] is defined by [L_nd]/[T_nd]
            >>> rvs_nd = np.array([[ 0.70239946, -0.43743181, -0.73375658,  0.15432556,  0.86144022,-0.35924967],[ 0.85918506, -0.5942174 , -0.89054218,  0.19227447,  0.98793658,-0.48574603]])
            >>> mu_nd = 1.5 
            >>> mee_nd = KeprvTrans.rv2mee(rvs_nd,mu_nd)
        Inputs:
            rvs -> [array-like,float] state vector
            mu -> [float] GM of the central attraction
            degrees -> [bool,optional,default=True] unit of L
        Outputs:
            mee -> [array-like] modified equinoctial elements in form of [p, f, g, h, k, L], where
                p -> [array-like,float] Semi-latus rectum, p = a*(1-e**2)
                f -> [array-like,float] x components of the eccentricity vector in the orbital frame, f = e * cos(Ω + ω)
                g -> [array-like,float] y components of the eccentricity vector in the orbital frame, g = e * sin(Ω + ω)
                h -> [array-like,float] x components of the node vector in the orbital frame, h = tan(i/2) * cos(Ω)
                k -> [array-like,float] y components of the node vector in the orbital frame, k = tan(i/2) * sin(Ω)
                L -> [array-like,float] True Longitude, [radians] or [deg], L = Ω + ω + ν
        """  
        mee = kep_rv_trans.rv2mee(rvs,mu,degrees) 

        return mee