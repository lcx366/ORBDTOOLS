from ..transform import orbele_trans                                                

class OrbeleTrans(object):
    """
    Class of transformation between classical orbital elements and non-singular orbital elements.
    """

    def nu_to_E(nu,ecc,degrees=True):
        """
        Transform to Eccentric Anomaly from True Anomaly.
    
        Usage:
            >>> from orbdtools import OrbeleTrans 
            >>> E = OrbeleTrans.nu_to_E(nu,ecc)
        Inputs:
            nu -> [array-like,float] True Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            E -> [array-like,float] Eccentric Anomaly, in (-π,π) or (-360,360)
        """
        E = orbele_trans.nu_to_E(nu,ecc,degrees)
        return E

    def nu_to_F(nu,ecc,degrees=True):
        """
        Transform to Hyperbolic Anomaly from True Anomaly.
    
        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> F = OrbeleTrans.nu_to_F(nu, ecc)
        Inputs:
            nu -> [array-like,float] True Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc > 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            F -> [array-like,float] Hyperbolic Anomaly, between -inf and inf, [radians] or [deg]
        """   
        F = orbele_trans.nu_to_F(nu,ecc,degrees) 
        return F

    def E_to_Me(E,ecc,degrees=True):
        """
        Transform to Mean Anomaly from Eccentric Anomaly.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> Me = OrbeleTrans.E_to_Me(E, ecc)
        Inputs:
            E -> [array-like,float] Eccentric Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            Me -> [array-like,float] Mean anomaly, [radians] or [deg]
        """
        Me = orbele_trans.E_to_Me(E,ecc,degrees)
        return Me

    def F_to_Mh(F,ecc,degrees=True):
        """
        Transform to Mean Anomaly from Hyperbolic Anomaly.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> Mh = OrbeleTrans.F_to_Mh(F, ecc)
        Inputs:
            F -> [array-like,float] Hyperbolic Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc > 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            Mh -> [array-like,float] Mean anomaly, [radians] or [deg]
        """  
        Mh = orbele_trans.F_to_Mh(F,ecc,degrees)
        return Mh

    def nu_to_Me(nu,ecc,degrees=True):
        """
        Transform to Mean anomaly from True Anomaly for ellipse trajectories.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> Me = nu_to_Me(nu, ecc)
        Inputs:
            nu -> [array-like,float] True Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            Me -> [array-like,float] Mean Anomaly, [radians] or [deg]
        """   
        Me = orbele_trans.nu_to_Me(nu,ecc,degrees) 
        return Me

    def nu_to_Mp(nu,degrees=True):
        """
        Transform to Mean anomaly from True Anomaly for parabolic trajectories.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> Mp = OrbeleTrans.nu_to_Mp(nu)
        Inputs:
            nu -> [array-like,float] True Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc = 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            Mp -> [array-like,float] Mean Anomaly, [radians] or [deg]
        """ 
        Mp = orbele_trans.nu_to_Mp(nu,degrees)
        return Mp

    def nu_to_Mh(nu,ecc,degrees=True):
        """
        Transform to Mean anomaly from True Anomaly for hyperbolic trahectories.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> Mh = OrbeleTrans.nu_to_Mh(nu, ecc)
        Inputs:
            nu -> [array-like,float] True Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc > 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            Mh -> [array-like,float] Mean Anomaly, [radians] or [deg]
        """ 
        Mh = orbele_trans.nu_to_Mh(nu,ecc,degrees)
        return Mh

    def Me_to_E(Me,ecc,degrees=True):
        """
        Transform to Eccentric Anomaly from Mean Anomaly.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> E = OrbeleTrans.Me_to_E(Me, ecc)
        Inputs:
            Me -> [array-like,float] Mean anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with 0 < ecc < 1 
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            E -> [array-like,float] Eccentric anomaly, [radians] or [deg]
        Notes:
            This uses a Newton iteration on the Kepler Equation.
        """ 
        E = orbele_trans.Me_to_E(Me,ecc,degrees)
        return E

    def Mh_to_F(Mh,ecc,degrees=True):
        """
        Transform to Hyperbolic Anomaly from Mean Anomaly.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> F = OrbeleTrans.Mh_to_F(Mh, ecc)
        Inputs:
            Mh -> [array-like,float] Mean anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc > 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            F -> [array-like,float] Hyperbolic anomaly
        Notes:
            This uses a Newton iteration on the Kepler Equation.
        """
        F = orbele_trans.Mh_to_F(Mh,ecc,degrees)
        return F

    def E_to_nu(E,ecc,degrees=True):
        """
        Transform to True anomaly from Eccentric Anomaly.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> nu = OrbeleTrans.E_to_nu(E, ecc)
        Inputs:
            E -> [array-like,float] Eccentric Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with 0 < ecc < 1 
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            nu -> [array-like,float] True Anomaly, in (-π,π) or (-360,360)
        """  
        nu = orbele_trans.E_to_nu(E,ecc,degrees) 
        return nu

    def F_to_nu(F, ecc,degrees=True):
        """
        Transform to True anomaly from Hyperbolic Anomaly.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> nu = OrbeleTrans.F_to_nu(F, ecc)
        Inputs:
            F -> [array-like,float] Hyperbolic Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc > 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            nu -> [array-like,float] True Anomaly, in (-π,π) or (-360,360)
        """ 
        nu = orbele_trans.F_to_nu(F, ecc,degrees)
        return nu

    def Me_to_nu(Me,ecc,degrees=True):
        """
        Transform to True Anomaly from Mean Anomaly for ellipse trajectories.
    
        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> nu = OrbeleTrans.Me_to_nu(Me, ecc)
        Inputs:
            Me -> [array-like,float] Mean Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with 0 < ecc < 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
        """
        nu = orbele_trans.Me_to_nu(Me,ecc,degrees)
        return nu

    def Mp_to_nu(Mp,degrees=True):
        """
        Transform to True Anomaly from Mean Anomaly for parabolic trajectories.
    
        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> nu = OrbeleTrans.Mp_to_nu(Mp)
        Inputs:
            Mp -> [array-like,float] Mean Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc = 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
        """   
        nu = orbele_trans.Mp_to_nu(Mp,degrees) 
        return nu

    def Mh_to_nu(Mh,ecc,degrees=True):
        """
        Transform to True Anomaly from Mean Anomaly for hyperbolic trajectories.
    
        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> nu = OrbeleTrans.Mh_to_nu(Mh, ecc)
        Inputs:
            Mh -> [array-like,float] Mean Anomaly, [radians] or [deg]
            ecc -> [array-like,float] Eccentricity with ecc > 1
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            nu -> [array-like,float] True Anomaly, in (-π,π) and (-360,360)
        """
        nu = orbele_trans.Mh_to_nu(Mh,ecc,degrees)
        return nu

    def coe2nse(a,ecc,inc,raan,argp,nu,degrees=True):
        """
        Convert to non-singular orbital elements from classical orbital elements for elliptic trajectories.
        The non-singular orbital elements exhibit no singularity for near-circular orbit, also known as the first kind of non-singular orbital elements.

        Usage: 
            >>> from orbdtools import OrbeleTrans
            >>> a, inc, raan, xi, eta, l = OrbeleTrans.coe2nse(a, ecc, inc, raan, argp, nu)
        Inputs:
            a -> [array-like,float] Semi-major axis
            ecc -> [array-like,float] Eccentricity
            inc -> [array-like,float] Inclination, [radians] or [deg]
            raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
            argp -> [array-like,float] Argument of perigee, [radians] or [deg]
            nu -> [array-like,float] True anomaly, [radians] or [deg]
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            a -> [array-like,float] Semi-major axis
            inc -> [array-like,float] Inclination, [radians] or [deg]
            raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
            xi -> [array-like,float] non-singular parameter, ξ = e * cosω
            eta -> [array-like,float] non-singular parameter, η = e * cosω
            l -> [array-like,float] non-singular parameter, [radians] or [deg], l = ω + M  
        """   
        nse = orbele_trans.coe2nse(a,ecc,inc,raan,argp,nu,degrees) 
        return nse

    def nse2coe(a,inc,raan,xi,eta,l,degrees=True):
        """
        Convert to classical orbital elements from non-singular orbital elements for elliptic trajectories.
        The non-singular orbital elements exhibit no singularity for near-circular orbit, also known as the first kind of non-singular orbital elements.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> a, ecc, inc, raan, argp, nu = OrbeleTrans.nse2coe(a, inc, raan, xi, eta, l)
        Inputs:
            a -> [array-like,float] Semi-major axis
            inc -> [array-like,float] Inclination, [radians] or [deg]
            raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
            xi -> [array-like,float] non-singular parameter, ξ = e * cosω
            eta -> [array-like,float] non-singular parameter, η = e * cosω
            l -> [array-like,float] non-singular parameter, [radians] or [deg], l = ω + M 
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            a -> [array-like,float] Semi-major axis
            ecc -> [array-like,float] Eccentricity
            inc -> [array-like,float] Inclination, [radians] or [deg]
            raan -> [array-like,float] Longitude of ascending node, [radians] or [deg]
            argp -> [array-like,float] Argument of perigee, [radians] or [deg]
            nu -> [array-like,float] True anomaly, [radians] or [deg]
        """   
        coe = orbele_trans.nse2coe(a,inc,raan,xi,eta,l,degrees) 
        return coe

    def coe2mee(a,ecc,inc,raan,argp,nu,degrees=True):
        """
        Convert to modified equinoctial orbital elements from classical orbital elements for elliptic trajectories. 
        The modified equinoctial orbital elements exhibit no singularity for near-circular orbit with inclinations close to 0 degrees. 
        It is also known as the second kind of non-singular orbital elements.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> p, f, g, h, k, L = OrbeleTrans.coe2mee(a, ecc, inc, raan, argp, nu)
        Inputs:
            a -> [array-like,float] Semi-major axis
            ecc -> [array-like,float] Eccentricity
            inc -> [array-like,float] Inclination, [radians] or [deg]
            omega -> [array-like,float] Longitude of ascending node, [radians] or [deg]
            argp -> [array-like,float] Argument of perigee, [radians] or [deg]
            nu -> [array-like,float] True anomaly, [radians] or [deg]
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            p -> [array-like,float] Semi-latus rectum, p = a*(1-e**2)
            f -> [array-like,float] Equinoctial parameter f, f = e * cos(Ω + ω)
            g -> [array-like,float] Equinoctial parameter g, g = e * sin(Ω + ω)
            h -> [array-like,float] Equinoctial parameter h, h = tan(i/2) * cos(Ω)
            k -> [array-like,float] Equinoctial parameter k, k = tan(i/2) * sin(Ω)
            L -> [array-like,float] Longitude, [radians] or [deg], L = Ω + ω + M
        """   
        mee = orbele_trans.coe2mee(a,ecc,inc,raan,argp,nu,degrees)  
        return mee

    def mee2coe(p,f,g,h,k,L,degrees=True):
        """
        Convert to classical orbital elements from modified equinoctial orbital elements for elliptic trajectories.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> a, ecc, inc, raan, argp, nu = OrbeleTrans.mee2coe(p, f, g, h, k, L)
        Inputs:
            p -> [array-like,float] Semi-latus rectum, p = a*(1-e**2)
            f -> [array-like,float] Equinoctial parameter f, f = e * cos(Ω + ω)
            g -> [array-like,float] Equinoctial parameter g, g = e * sin(Ω + ω)
            h -> [array-like,float] Equinoctial parameter h, h = tan(i/2) * cos(Ω)
            k -> [array-like,float] Equinoctial parameter k, k = tan(i/2) * sin(Ω)
            L -> [array-like,float] Longitude, [radians] or [deg], L = Ω + ω + M
            degrees -> [bool,optional,default=True] unit for angular variables 
        Outputs:
            a -> [array-like,float] Semi-major axis
            ecc -> [array-like,float] Eccentricity
            inc -> [array-like,float] Inclination, [radians] or [deg]
            omega -> [array-like,float] Longitude of ascending node, [radians] or [deg]
            argp -> [array-like,float] Argument of perigee, [radians] or [deg]
            nu -> [array-like,float] True anomaly, [radians] or [deg]
        """
        coe = orbele_trans.mee2coe(p,f,g,h,k,L,degrees)
        return coe

    def mean2osculating(mean_ele,epoch,meanref='TEME',oscuref='TEME',degrees=True):
        """
        Convert mean orbital elements to osculating orbital elements using sgp4/sdp4.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> from astropy.time import Time
            >>> mean_ele = [7000,0.01,50,100,30,210] # in form of [a, e, i, Ω, ω, v]
            >>> epoch = Time('2022-06-07T08:09:12.345')
            >>> oscu_ele = OrbeleTrans.mean2osculating(mean_ele,epoch)
        Inputs:
            mean_ele -> [list or array of float] mean elements for sgp4/sdp4
            epoch -> Object of class Astropy Time
            meanref -> [str,optional,default='TEME'] reference frame bound by the mean elements
            oscuref -> [str,optional,default='TEME'] reference frame bound by the osculating elements
            degrees -> [bool,optional,default=True] unit of the angular variable of orbital elements
        Outputs:
            oscu_ele -> [list or array of float] osculating elements for sgp4/sdp4
        """    
        oscu_ele = orbele_trans.mean2osculating(mean_ele,epoch,meanref,oscuref,degrees)
        return oscu_ele

    def osculating2mean(oscu_ele,epoch,oscuref='TEME',meanref='TEME',degrees=True):
        """
        Convert osculating orbital elements to mean orbital elements using sgp4/sdp4.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> from astropy.time import Time
            >>> oscu_ele = [7000,0.01,50,100,30,210]
            >>> epoch = Time('2022-06-07T08:09:12.345')
            >>> mean_ele = OrbeleTrans.mean2osculating(oscu_ele,epoch)
        Inputs:
            oscu_ele -> [list or array of float] osculating elements for sgp4/sdp4
            epoch -> Object of class Astropy Time
            oscuref -> [str,optional,default='TEME'] reference frame bound by the osculating elements
            meanref -> [str,optional,default='TEME'] reference frame bound by the mean elements
            degrees -> [bool,optional,default=True] unit of the angular variable of orbital elements
        Outputs:
            mean_ele -> [list or array of float] mean elements for sgp4/sdp4
        """ 
        mean_ele = orbele_trans.osculating2mean(oscu_ele,epoch,oscuref,meanref,degrees) 
        return mean_ele

    def coe_trans(trans_matrix,coe_from,degrees=True):
        """
        Converting classical orbital elements between two reference frames, especially 'TEME' and 'ICRF'. The reference frame must be defined as hand-right.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> coe_from = [6.9974e+03,1.0673e-02,4.9991e+01,1.0002e+02,3.3132e+01,2.0688e+02] # in TEME
            >>> epoch = Time('2022-06-07T08:09:12.345')
            >>> trans_matrix = FrameTrans.gcrf_teme_mat(epoch)
            >>> teme2gcrf_mat = trans_matrix.teme2gcrf_mat # transformation matrix from TEME to GCRF
            >>> coe_to = OrbeleTrans.coe_trans(teme2gcrf_mat,coe_from)
            >>> print(coe_to)
        Inputs:
            trans_matrix -> [nx3x3 array-like] single or mutiple transformation matrix
            coe_from -> [nx6] classical orbital elements in source reference frame
            degrees -> [bool,optional,default=True] unit of the angular variable of orbital elements
        Outputs:
            coe_to -> [nx6] classical orbital elements in target reference frame  
        """ 
        coe_to = orbele_trans.coe_trans(trans_matrix,coe_from,degrees)

        return coe_to

    def rv_trans(trans_matrix,rv_from):
        """
        Converting orbital state vector between two reference frames.

        Usage:
            >>> from orbdtools import OrbeleTrans
            >>> from orbdtools import FrameTrans
            >>> from astropy.time import Time
            >>> epoch = Time('2022-06-07T08:09:12.345')
            >>> trans_matrix = FrameTrans.gcrf_teme_mat(epoch)
            >>> teme2gcrf_mat = trans_matrix.teme2gcrf_mat # transformation matrix from TEME to GCRF 
            >>> rv_from = np.array([4.4836e+03,-2.7941e+03,-4.6840e+03,1.2189,6.8128,-2.8404]) # in unit of [km,km/s] in TEME
            >>> rv_to = OrbeleTrans.rv_trans(teme2gcrf_mat,rv_from)
            >>> print(rv_to)

            >>> rv_to = OrbeleTrans.rv_trans(trans_matrix,rv_from)
        Inputs:
            trans_matrix -> [nx3x3 array-like] single or mutiple transformation matrix
            rv_from -> [nx6] orbital state vector in source reference frame
        Outputs:
            rv_to -> [nx6] orbital state vector in target reference frame  
        """
        rv_to = orbele_trans.rv_trans(trans_matrix,rv_from)
        return rv_to