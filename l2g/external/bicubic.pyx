#distutils: language = c++
#cython: language_level = 3

import numpy as np
cimport numpy as np

# External CPP class
from l2g.external.bicubic cimport BICUBIC_INTERP, PyBicubic
from libcpp.vector cimport vector


cdef class PyBicubic:
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
    def __cinit__(self):
        self.c_bicubic = new BICUBIC_INTERP()
        # self.c_BI_DATA = BI_DATA()

    def __init__(self, x, y, f):
        """Initialize the class. If we provide the three arrays then
        """
        self.prepared: bool = False

        if x is None or y is None or f is None:
            return
        else:
            self.setArrays(x, y, f)

    def setArrays(self, x: np.ndarray, y: np.ndarray, f: np.ndarray):
        if not isinstance(x, np.ndarray):
            if isinstance(x, list):
                x = np.array(x)
            else:
                raise Exception("Argument x is not a numpy array")

        if not isinstance(y, np.ndarray):
            if isinstance(y, list):
                y = np.array(x)
            else:
                raise Exception("Argument y is not a numpy array")

        if not isinstance(f, np.ndarray):
            if isinstance(f, list):
                f = np.array(f)
            else:
                raise Exception("Argument f is not a numpy array")

        if x.ndim != 1:
            raise Exception("Argument x is not a 1D array")

        if y.ndim != 1:
            raise Exception("Argument y is not a 1D array")

        if f.ndim != 2:
            raise Exception("Argument f is not a 2D array")

        if x.shape[0] != f.shape[1]:
            raise Exception(f"Shape mismatch with x ({x.shape[0]}) and f ({f.shape[1]})")

        if y.shape[0] != f.shape[0]:
            raise Exception(f"Shape mismatch with y ({y.shape[0]}) and f ({f.shape[0]})")

        cdef:
            vector[double] vx, vy
            vector[vector[double]] vf
            vector[double] buff

        vx.clear()
        vy.clear()
        vf.clear()

        for i in range(x.size):
            vx.push_back(x[i])

        for i in range(y.size):
            vy.push_back(y[i])

        # Reshape the array
        for i in range(f.shape[0]):
            buff.clear()
            for j in range(f.shape[1]):
                buff.push_back(f[i][j])
            vf.push_back(buff)

        self.c_bicubic.setArrays(vx, vy, vf)
        self.c_bicubic.populateContext(&self.c_BI_DATA)

    def __call__(self, x, y):
        """Passed x,y can be an array.
        """
        nx = 1
        if isinstance(x, list):
            nx = len(x)

        ny = 1
        if isinstance(y, list):
            ny = len(y)

        if nx == ny:
            if nx == 1:
                return self.getValues(x, y)
            else:
                out_val = []
                out_valdx  = []
                out_valdy = []
                for i in range(nx):
                    val, valdx, valdy = self.getValues(x[i], y[i])
                    out_val.append(val)
                    out_valdx.append(valdx)
                    out_valdy.append(valdy)
                return out_val, out_valdx, out_valdy
        else:
            raise Exception(f"Mismatch in number of elements in x and y: {nx} {ny}")

    def isPrepared(self):
        return self.prepared

    def getValues(self, double x, double y):
        """Returns the interpolated value of the function at point (x, y)

        Arguments:
            x: X component of the point
            y: Y component of the point

        Returns:
            val: Value of the function at point (x,y)
            valdx: Value of the first derivative (df/dx) at (x, y)
            valdy: Value of the first derivative (df/dy) at (x, y)
        """

        self.c_BI_DATA.r = x
        self.c_BI_DATA.z = y
        self.c_bicubic.getValues(&self.c_BI_DATA)

        return self.c_BI_DATA.val, self.c_BI_DATA.valdx, self.c_BI_DATA.valdy

    def getSecondDerivativeValues(self, double x, double y):
        """Returns the interpolated value of the function at point (x, y)

        Arguments:
            x: X component of the point
            y: Y component of the point

        Returns:
            valdxdy: First mixed derivative value of the function at point (x,y)
            valdxdx: Value of the second derivative (d^2f/dx^2) at (x, y)
            valdydy: Value of the second derivative (d^2f/dy^2) at (x, y)
        """
        self.c_BI_DATA.r = x
        self.c_BI_DATA.z = y
        self.c_bicubic.getSecondDerivativeValues(&self.c_BI_DATA)

        return self.c_BI_DATA.val, self.c_BI_DATA.valdx, self.c_BI_DATA.valdy
