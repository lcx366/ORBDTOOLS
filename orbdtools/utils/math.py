import numpy as np

def Matrix_dot_Vector(matrix,vector):
    """
    Compute the dot product of matrices and a vectors.

    Usage:
        >>> import numpy as np
        >>> matrix = np.arange(1296).reshape(8,18,3,3)
        >>> vector = np.arange(432).reshape(8,18,3)
        >>> vector_trans = Matrix_dot_Vector(matrix,vector)
    Inputs:
        matrix -> [array like] Multi-dimensional Matrix
        vector -> [array like] Multi-dimensional Vector
    Outputs:
        vector_trans -> [array like] Transformed vectors          
    """
    matrix,vector = np.array(matrix),np.array(vector)
    if matrix.ndim == 2 and vector.ndim == 1:
        vector_trans = matrix @ vector
    elif matrix.ndim == 3 and vector.ndim == 2:
        vector_trans = np.squeeze(matrix @ vector[:,:,None])
    elif matrix.ndim == 4 and vector.ndim == 3:	
        vector_trans = np.squeeze(matrix @ vector[:,:,:,None])
    else:
        raise Exception('The dimensions of the matrix and vector do not match.')
        
    return vector_trans

def getangle(sin,cos):
    """
    Calculate the angle from values of sin and cos
    """
    return np.arctan2(sin,cos)    