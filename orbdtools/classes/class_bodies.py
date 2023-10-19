import numpy as np

class Body(object):
    
    def __init__(self,info):
        """
        Initialize a celestial body.
        """
        L_nd_unit,mu_nd_unit = info['Re'],info['mu'] # Non-dimensional length unit [L_nd] and GM unit [mu_nd]
        T_nd_unit = np.sqrt(L_nd_unit**3/mu_nd_unit) # Non-dimensional time unit [T_nd]
        v_nd_unit = L_nd_unit/T_nd_unit
        dict_keys = ['_L','_mu','_T','_v']
        dict_values = [L_nd_unit,mu_nd_unit,T_nd_unit,v_nd_unit]
        add_info = dict(zip(dict_keys, dict_values))
        info.update(add_info)

        for key in info.keys():
            setattr(self, key, info[key])
            
    def __repr__(self):
        """
        Returns a more information-rich string representation of the Body object.
        """
        return f"<Body object: {self.name}({self.symbol})>" 

    def create_body(info):
        """
        Create a celestial body.

        Usage:
            >>> from orbdtools import Body
            >>> info_earth = {'name':"Earth",'symbol':"\u2641",
            >>> 'mu':398600.4418, # GM, [km^3/s^2]
            >>> 'Re':6378.137, # Equatorial radius of the Earth from WGS84 ellipsoid model, [km]
            >>> 'f':1/298.257223563, # flattening
            >>> 'R':6371.0008, # volumetric radius, [km] 
            >>> 'J2':1.08262982e-3}
            >>> earth = Body(info_earth)
        Inputs:
            info -> [dictionary] Necessary input information used to create a celestial body
        Outputs:
            body -> [object of class Body] Celestial body created   
        """
        body = Body(info) 
        return body   

    def from_name(bodyname):  
        """
        Load a celestial body.

        Usage:
            >>> from orbdtools import Body
            >>> earth = Body.from_name('Earth')
        Inputs:
            bodyname -> [str] Name of the celestial body to load. Currently, available options include 'Earth', 'Moon', 'Sun', 'Venus', 'Jupiter', and 'Mars'.
        Outputs:
            body -> [object of class Body] Celestial body loaded
        """
        if bodyname == 'Earth':  
            info = {'name':"Earth",
              'symbol':"\u2641",
              'mu':398600.4418, # GM, [km^3/s^2]
              'Re':6378.137, # Equatorial radius of the Earth from WGS84 ellipsoid model, [km]
              'f':1/298.257223563, # flattening
              'R':6371.0008, # volumetric radius, [km] 
              'J2':1.08262982e-3,
              'n_sunsync':1.99096871e-7, # the precession rate of Sun-synchronous orbit in rad/s
              'stellar_year':365.2564 # mean solar day
              }
        elif bodyname == 'Moon':          
            info = {'name':"Moon",
             'symbol':"\u263E",
             'mu':4902.8001,   
             'Re':1738.1,
             'R':1737.4,
             'f':0.0012,
             'J2':2.0335425e-4
              }
        elif bodyname == 'Sun':
            info = {'name':"Sun",
             'symbol':"\u2609",
             'mu':1.32712440018e11,
             'Re':695700,
             'R':695694,
             'f':9e-6,
             'J2':2.211e-7
              }
        elif bodyname == 'Venus':      
            info = {'name':"Venus",
             'symbol':"\u2640",
             'mu':324859,   
             'R':6051.8,
             'f':0,
             'J2':4.4044e-6
              }
        elif bodyname == 'Jupiter':
            info = {'name':"Jupiter",
             'symbol':"\u2643",
             'mu':1.26686534e8,
             'Re':71492,
             'R':69911.3,
             'f':0.064874,
             'J2':0.01475
              }
        elif bodyname == 'Mars':      
            info = {'name':"Mars",
             'symbol':"\u2642",
             'mu':42828.37,
             'Re':3396.19,
             'R':3389.5,
             'f':0.00589,
             'J2':0.0019555
              }
        else:
            raise Exception("{:s} is not in the depository, and you may try to CREATE it by 'Body.create_body()'.".format(bodyname))  
        return Body(info)                         