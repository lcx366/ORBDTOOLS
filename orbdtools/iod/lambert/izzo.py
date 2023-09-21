import numpy as np
from numpy.linalg import norm

from ...transform.kep_rv_trans import rv2coe

def _find_xy(ll, T, M, numiter, rtol):
    """
    Computes all x, y for given number of revolutions.
    """
    # For abs(ll) == 1 the derivative is not continuous

    M_max = np.floor(T / np.pi)
    T_00 = np.arccos(ll) + ll * np.sqrt(1 - ll ** 2)  # T_xM

    # Refine maximum number of revolutions if necessary
    if T < T_00 + M_max * np.pi and M_max > 0:
        _, T_min = _compute_T_min(ll, M_max, numiter, rtol)
        if T < T_min: M_max -= 1

    # Check if a feasible solution exist for the given number of revolutions
    # This departs from the original paper in that we do not compute all solutions
    if M > M_max: raise ValueError("No feasible solution, try lower M")

    # Initial guess
    for x_0 in _initial_guess(T, ll, M):
        # Start Householder iterations from x_0 and find x, y
        x = _householder(x_0, T, ll, M, rtol, numiter)
        y = _compute_y(x, ll)
        yield x, y # generate an iterators

def _compute_T_min(ll, M, numiter, rtol):
    """
    Compute minimum T.
    """
    if ll == 1:
        x_T_min = 0
        T_min = _tof_equation(x_T_min, 0, ll, M)
    else:
        if M == 0:
            x_T_min = np.inf
            T_min = 0.0
        else:
            # Set x_i > 0 to avoid problems at ll = -1
            x_i = 0.1
            T_i = _tof_equation(x_i, 0, ll, M)
            x_T_min = _halley(x_i, T_i, ll, rtol, numiter)
            T_min = _tof_equation(x_T_min, 0, ll, M)

    return x_T_min, T_min   

def _tof_equation(x, T0, ll, M):
    """
    Time of flight equation.
    """
    y = _compute_y(x, ll)
    return _tof_equation_y(x, y, T0, ll, M)     

def _halley(p0, T0, ll, tol, maxiter):
    """
    Find a minimum of time of flight equation using the Halley method.
    Notes:
    This function is private because it assumes a calling convention specific to
    this module and is not really reusable.
    """
    flag = False
    for j in range(maxiter):
        y = _compute_y(p0, ll)
        fder = _tof_equation_p(p0, y, T0, ll)
        fder2 = _tof_equation_p2(p0, y, T0, fder, ll)
        if fder2 == 0: raise Exception("Derivative was zero")
        fder3 = _tof_equation_p3(p0, y, T0, fder, fder2, ll)

        # Halley step (cubic)
        p = p0 - 2 * fder * fder2 / (2 * fder2 ** 2 - fder * fder3)

        if np.abs(p - p0) < tol: 
            flag = True
            break
        else:    
            p0 = p

    if flag:  
       return p   
    else:      
        raise RuntimeError("Failed to converge")    

def _tof_equation_y(x, y, T0, ll, M):
    """
    Time of flight equation with externally computated y.
    """
    if M == 0 and np.sqrt(0.6) < x < np.sqrt(1.4):
        eta = y - ll * x
        S_1 = (1 - ll - x * eta) * 0.5
        Q = 4 / 3 * hyp2f1b(S_1)
        T_ = (eta ** 3 * Q + 4 * ll * eta) * 0.5
    else:
        psi = _compute_psi(x, y, ll)
        T_ = np.divide(np.divide(psi + M * np.pi, np.sqrt(np.abs(1 - x ** 2))) - x + ll * y,(1 - x ** 2))

    return T_ - T0

def _tof_equation_p(x, y, T, ll):
    # TODO: What about derivatives when x approaches 1?
    res = (3 * T * x - 2 + 2 * ll ** 3 * x / y) / (1 - x ** 2)  
    return res
    
def _tof_equation_p2(x, y, T, dT, ll):
    res = (3 * T + 5 * x * dT + 2 * (1 - ll ** 2) * ll ** 3 / y ** 3) / (1 - x ** 2)
    return res


def _tof_equation_p3(x, y, _, dT, ddT, ll):
    res = (7 * x * ddT + 8 * dT - 6 * (1 - ll ** 2) * ll ** 5 * x / y ** 5) / (1 - x ** 2)       
    return res

def _compute_y(x, ll):
    """
    Computes y.
    """
    y = np.sqrt(1 - ll ** 2 * (1 - x ** 2))
    return y      

def _compute_psi(x, y, ll):
    """
    Computes psi.
    The auxiliary angle psi is computed using Eq.(17) by the appropriate
    inverse function
    """
    if -1 <= x < 1:
        # Elliptic motion
        # Use arc cosine to avoid numerical errors
        return np.arccos(x * y + ll * (1 - x ** 2))
    elif x > 1:
        # Hyperbolic motion
        # The hyperbolic sine is bijective
        return np.arcsinh((y - x * ll) * np.sqrt(x ** 2 - 1))
    else:
        # Parabolic motion
        return 0   

def hyp2f1b(x):
    """
    Hypergeometric function 2F1(3, 1, 5/2, x), see [Battin].

    Notes:
    More information about hypergeometric function can be checked at
    https://en.wikipedia.org/wiki/Hypergeometric_function
    """
    if x >= 1:
        return np.inf
    else:
        res = term = 1
        i = 0
        while term > 2e-16:
            term *= (3 + i) * (1 + i) / (5 / 2 + i) * x / (i + 1)
            res += term
            i += 1 
        return res    

def _initial_guess(T, ll, M):
    """
    Initial guess.
    """
    if M == 0:
        # Single revolution
        T_0 = np.arccos(ll) + ll * np.sqrt(1 - ll ** 2) + M * np.pi  # Equation 19
        T_1 = 2 * (1 - ll ** 3) / 3  # Equation 21
        if T >= T_0:
            x_0 = (T_0 / T) ** (2 / 3) - 1
        elif T < T_1:
            x_0 = 5 / 2 * T_1 / T * (T_1 - T) / (1 - ll ** 5) + 1
        else:
            # This is the real condition, which is not exactly equivalent
            # elif T_1 < T < T_0
            # piecewise equation right after expression (30) in the original paper is incorrect
            # See https://github.com/poliastro/poliastro/issues/1362
            x_0 = 2**(np.log(T / T_0) / np.log(T_1 / T_0)) -1

        return [x_0]
    else:
        # Multiple revolution
        x_0l = (((M * np.pi + np.pi) / (8 * T)) ** (2 / 3) - 1) / (((M * np.pi + np.pi) / (8 * T)) ** (2 / 3) + 1)
        x_0r = (((8 * T) / (M * np.pi)) ** (2 / 3) - 1) / (((8 * T) / (M * np.pi)) ** (2 / 3) + 1)

        return [x_0l, x_0r]  
        
def _householder(p0, T0, ll, M, tol, maxiter):
    """
    Find a zero of time of flight equation using the Householder method.
    Notes:
    This function is private because it assumes a calling convention specific to
    this module and is not really reusable.
    """
    flag = False
    for j in range(maxiter):
        y = _compute_y(p0, ll)
        fval = _tof_equation_y(p0, y, T0, ll, M)
        T = fval + T0
        fder = _tof_equation_p(p0, y, T, ll)
        fder2 = _tof_equation_p2(p0, y, T, fder, ll)
        fder3 = _tof_equation_p3(p0, y, T, fder, fder2, ll)

        # Householder step (quartic)
        p = p0 - fval * ((fder ** 2 - fval * fder2 / 2)/ (fder * (fder ** 2 - fval * fder2) + fder3 * fval ** 2 / 6))

        if np.abs(p - p0) < tol:
            flag = True
            break
        else:
            p0 = p

    if flag:
        return p  
    else:          
        raise RuntimeError("Failed to converge")   
    
def _reconstruct(x, y, r1, r2, ll, gamma, rho, sigma):
    """
    Reconstruct solution velocity vectors.
    """
    V_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1
    V_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2
    V_t1 = gamma * sigma * (y + ll * x) / r1
    V_t2 = gamma * sigma * (y + ll * x) / r2
    return V_r1, V_r2, V_t1, V_t2                  

def izzo_iod(mu,r1_vec,r2_vec,tof,tm=1,M=0,degrees=True,numiter=35,rtol=1e-8):
    """
    Solve the Lambert's problem using the Izzo algorithm.

    Inputs:
        mu -> [float] GM of the center of attraction
        r1_vec -> Initial position in normalized Length unit
        r2_vec -> Final position in normalized Length unit
        tof -> Time of flight in normalized Time unit
    Parameters:    
        M -> [int,default=0] Number of full revolutions
        tm -> [int, default=1] Transfer mode, tm = 1 for short way and -1 for long way.
        numiter -> [int,default=35] Maximum number of iterations
        rtol -> [float, default=1e-8] Relative tolerance of the algorithm
    Outputs:
        v1, v2 -> [tuple of float] Pair of velocity solutions
    """
    eles = []

    # Check collinearity of r1 and r2
    if not np.cross(r1_vec, r2_vec).any():
        raise ValueError("Lambert solution cannot be computed for collinear vectors")

    # Chord
    c_vec = r2_vec - r1_vec
    c, r1, r2 = norm([c_vec,r1_vec,r2_vec],axis=1)

    # Semiperimeter
    s = (r1 + r2 + c)/2

    # unit versors
    r1_uec, r2_uec = r1_vec / r1, r2_vec / r2

    # For shot way 
    h_uec = np.cross(r1_uec, r2_uec)
    ll = np.sqrt(1 - min(1, c / s))

    # for long way
    if tm == -1:
        h_uec = -h_uec
        ll = -ll
    h_uec = h_uec / norm(h_uec)  # Fixed from paper

    t1_uec, t2_uec = np.cross(h_uec, r1_uec), np.cross(h_uec, r2_uec)  # Fixed from paper

    # Non dimensional time of flight
    T = np.sqrt(2 * mu / s ** 3) * tof

    # Find solutions
    xy = _find_xy(ll, T, M, numiter, rtol)

    # Reconstruct
    gamma = np.sqrt(mu * s / 2)
    rho = (r1 - r2) / c
    sigma = np.sqrt(1 - rho ** 2)

    for x,y in xy:
        # Calculate velocity projected on radial and transverse(along-track) direction
        vr1, vr2, vt1, vt2 = _reconstruct(x, y, r1, r2, ll, gamma, rho, sigma)
        v1_vec = vr1 * r1_uec + vt1 * t1_uec
        v2_vec = vr2 * r2_uec + vt2 * t2_uec
        ele = rv2coe(np.hstack([r1_vec,v1_vec]),mu,degrees) 
        eles.append(ele)
    eles = np.array(eles)    
    return eles