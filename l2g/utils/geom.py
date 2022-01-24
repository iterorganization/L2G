import numpy as np

def getBaryCenter(points, coord='cyl', mul=1e-3):
    """Get barycenter of triangle.

    Arguments:
        points (np.float): 1D array of floats of size 9
        coord (str): Type of coordinates to return
        mul (float): Multiplier for dimension. For conversion from meter to mm,
                     or vice versa
    """

    bc = np.empty(3, dtype=np.float)

    bc[0] = np.sum(points[0::3])
    bc[1] = np.sum(points[1::3])
    bc[2] = np.sum(points[2::3])
    bc /= 3


    if coord == 'cyl':
        phi = np.arctan2(bc[1], bc[0])
        r = np.sqrt(bc[0]**2 + bc[1]**2)
        return mul * r, mul * bc[2], phi
    else:
        return mul * bc

def getTriangleNormal(points):
    """Get normal of triangle. Orientation of the points must be positive!

    Therefor the points list must be:
        [p1[0], p1[1], p1[2], p2[0]....]

    Arguments:
        points (np.float): 1D array of floats of size 9
    """
    normalVec = np.cross(points[3:6] - points[:3], points[6:] - points[:3])
    normalVec /= np.linalg.norm(normalVec)
    return normalVec

def angleBetweenVectors(v1, v2):
    """Return the angle between two vectors
    """
    dot = np.dot(v1, v2)
    sq1 = np.dot(v1, v1)
    sq2 = np.dot(v2, v2)

    angle = np.arccos(dot / np.sqrt(sq1 * sq2))
    return angle
