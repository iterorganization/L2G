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

def rotatePointsAroundAxis(points: np.ndarray, p1: np.ndarray, p2: np.ndarray,
    theta: float):
    """Rotate the argument points around the axis |p2 - p1| for the angle
    theta.

    Arguments:
        points (np.ndarray): Array of points with the shape ((N_points, 3)).
        p1 (np.ndarray): Start point of the axis
        p2 (np.ndarray): End point of the axis
        theta (float): angle in degrees.

    Returns:
        rotated_points (np.ndarray)
    """

    # Copy the points and translate it to p1.
    out_points = np.zeros(points.shape)
    out_points[:] = points - p1

    # Calculate the normed vector
    diff_p = p2 - p1
    u = diff_p / np.linalg.norm(diff_p)

    # Construct the rotational matrix
    R = np.zeros((3, 3))

    cos_th = np.cos(np.deg2rad(theta))
    sin_th = np.sin(np.deg2rad(theta))
    moin_cos_th = 1 - cos_th
    ux = u[0]
    uy = u[1]
    uz = u[2]
    R[0, 0] = cos_th + ux * ux * moin_cos_th
    R[1, 0] = uy * ux * moin_cos_th + uz * sin_th
    R[2, 0] = uz * ux * moin_cos_th - uy * sin_th

    R[0, 1] = ux * uy * moin_cos_th - uz * sin_th
    R[1, 1] = cos_th + uy * uy * moin_cos_th
    R[2, 1] = uz * uy * moin_cos_th + ux * sin_th

    R[0, 2] = ux * uz * moin_cos_th + uy * sin_th
    R[1, 2] = uy * uz * moin_cos_th - ux * sin_th
    R[2, 2] = cos_th + uz * uz * moin_cos_th

    # Rotate the points
    # out_points = R@out_points.T
    out = (R@out_points.T).T

    # We need to fix the stride. Best way is to copy all values...
    # No idea otherwise how to fix the array stride.
    out_points[:] = out
    # out_points = out_points.T

    # Re-translate it to correct place.
    # out_points += p1
    return out_points + p1
