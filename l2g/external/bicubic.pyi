import numpy as np

class PyBicubic:
    """Function that wraps an implementation of the Bicubic interpolation
    algorithm for 2D functions defined on a rectilinear grid.


    .. code-block:: python

       import l2g.external.bicubic

       # Setup input data arrays
       x = ... # 1D array
       y = ... # 1D array
       f = ... # 2D array

       obj = l2g.external.bicubic.PyBicubic(x, y, f)

       # Query for values
       xp: float
       yp: float

       v, vdx, vdy = obj.getValues(xp, yp)
       v, vdx, vdy = obj(xp, yp)
       av, avdx, avdx = obj([xp1, xp2, ...], [yp1, yp2, ...])

       # Second derivative values
       vdxdy,vdxdx, vdydy = obj.getSecondaryDerivatives(xp, yp)

    """
    def __init__(self, x: np.ndarray, y: np.ndarray, f: np.ndarray) -> None: ...
    def setArrays(self, x: np.ndarray, y: np.ndarray, f: np.ndarray) -> None: ...
    def isPrepared(self) -> bool: ...
    def __call__(self, x: float, y: float) -> tuple[float, float, float]: ...
    def getValues(self, x: float, y: float) -> tuple[float, float, float]:
        """Returns the interpolated value of the function at point (x, y)

        Arguments:
            x: X component of the point
            y: Y component of the point

        Returns:
            val: Value of the function at point (x,y)
            valdx: Value of the first derivative (df/dx) at (x, y)
            valdy: Value of the first derivative (df/dy) at (x, y)
        """
    def getSecondDerivativeValues(self, x: float, y: float) -> tuple[float, float, float]:
        """Returns the interpolated value of the function at point (x, y)

        Arguments:
            x: X component of the point
            y: Y component of the point

        Returns:
            valdxdy: First mixed derivative value of the function at point (x,y)
            valdxdx: Value of the second derivative (d^2f/dx^2) at (x, y)
            valdydy: Value of the second derivative (d^2f/dy^2) at (x, y)
        """