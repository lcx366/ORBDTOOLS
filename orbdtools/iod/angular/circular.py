from ..radar.gibbs import gibbs,herrick_gibbs,gibbs_assist
from ..common import slant,get_a0

def solve_circular(mu,tof,tau,los3p,xyz_site3p,degrees=True):
    """
    Estimate the classical orbital elements at Median epoch from optical angle-only measurements using Near-Circular Orbit Hypothesis method.

    Usage:
        >>> ele = solve_circular(mu,tof,tau,los3p,xyz_site3p)
    Inputs:
        mu -> [float] GM of the central body of attraction
        tof -> [float] Time of flight
        tau -> [tuple of float] (t2-t1,t3-t2), where (t1,t2,t3) represents the start, midpoint, and end times of the observation
        los3p -> [2D array with shape of 3x3] Line-Of-Sight(LOS) vector of the space object relative to the site at (t1,t2,t3)
        xyz_site3p -> [2D array with shape of nx3] Cartesian coordinates of the site at (t1,t2,t3)
        degrees -> [bool,optional,default=True] Unit of angular variables in classical orbital elements estimated. If True, angular variables are in [deg], otherwise in [rad].
    Outputs:
        ele -> [list] Classical orbital elements with values as follows
            a -> [float] Semi-major axis
            ecc -> [float] Eccentricity
            inc -> [float] Inclination, [rad] or [deg]
            raan -> [float] Longitude of ascending node, [rad] or [deg]
            argp -> [float] Argument of perigee, [rad] or [deg]
            nu -> [float] True anomaly, [rad] or [deg]
    references: 
        张晓祥,吴连大,熊建宁.空间目标的圆轨道跟踪法[J].天文学报, 2003, 44(4):11.DOI:10.3321/j.issn:0001-5245.2003.04.010.    
    """
    # Calculate the semi-major axis of the circular orbit
    a = get_a0(mu,los3p,xyz_site3p,tof)

    # compute the cartesian coordinates of the space object in GCRF at (t1,t2,t3)
    rho,r_vec = slant(a,los3p,xyz_site3p)

    # Judge whether to use the method Gibbs or Herrick-Gibbs to determine the initial orbit based on three position vectors.
    method,mode,nu21,nu32,nu31 = gibbs_assist(r_vec,degrees)
    
    if method == 'Gibbs':
        ele = gibbs(mu,r_vec,degrees)
    elif method == 'Herrick-Gibbs':
        ele = herrick_gibbs(mu,r_vec,tof,tau,degrees) 
    
    return ele