# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np

from libc.math cimport atan2
from libc.stdio cimport printf

from l2g.external.bicubic cimport PyBicubic
from l2g.external.rkf45 cimport PyRKF45FLT

cdef class PyRKF45FLT:
    """Class that wraps a C++ implementation of the RKF4(5) method for solving
    the fieldline tracing equations in an axisymmetrical Cylindrical coordinate
    system.

    .. note::

       This is not an implementation that can solve a general system of PDEs.
       This class is meant for easily plotting specific fieldlines.

    In other words it is a class that follows contours of a 2D function,
    defined in a (R, Z) Cylindrical space. Based on the

    .. code_block:: python

       from l2g.external.bicubic import PyBicubic
       from l2g.external.rkf45 import PyRKF45FLT


       int_obj = PyBicubic(x, y, f)

       rkf45_obj = PyRKF45FLT()

       # Set the contour data
       rkf45_obj.setInterpolator(int_obj)
       # Essentially set the constant value F = B_t * R, from which B_t is
       # evaluated as F / R.
       rkf45_obj.set_vacuum_fpol(vacuum_fpol)


       # Start following
       r: float
       z: float
       th: float
       time_step: float
       nr: float
       nz: float
       nth: float

       while 1:
           nr, nz, nth = rkf45_obj.run_step(r, z, th, time_step)
    """

    def __cinit__(self):
        self.c_rkf45 = new RKF45()
        self.c_rkf45.set_omp_thread(0)

    def __init__(self):
        self.c_relerr = 1e-4
        self.c_abserr = 1e-4

    def setInterpolator(self, PyBicubic obj):
        """Sets the interpolator object that contains the data for following
        contours.

        Arguments:
            obj (PyBicubic): Bicubic interpolator used for obtaining the
                gradient values of the interpolated data.
        """
        self.c_rkf45.set_interpolator(obj.c_bicubic)

    def set_r_rmove(self, double r_move):
        """Sets the radial shift of the data.

        Arguments:
            r_move (float): Radial shift.
        """
        self.c_rkf45.set_r_move(r_move)

    def set_z_rmove(self, double z_move):
        """Sets the vertical shift of the data.

        Arguments:
            r_move (float): Vertical shift.
        """
        self.c_rkf45.set_z_move(z_move)

    def set_vacuum_fpol(self, double vacuum_fpol):
        """Sets the constant F = B_t R, used for calculating the toroidal
        component of the magnetic field as B_t = F / R.

        Arguments:
            vacuum_fpol (float): Poloidal current value.
        """
        self.c_rkf45.set_vacuum_fpol(vacuum_fpol)

    def run_step(self, double r, double z, double th, double time_step):
        """Perform one step of RKF4(5) from point (r, z, th) in the Cylindrical
        coordinate system.

        Arguments:
            r (float): Radial position
            z (float): Vertical position
            th (float): Toroidal angle
            time_step (float): Size of time_step

        Returns:
            r (float): New radial position
            z (float): Vertical position
            th (float): New toroidal angle
        """

        cdef:
            double y[2]
            double yp[2]
            double relerr, abserr, time, new_time
            int flag
        relerr = self.c_relerr
        abserr = self.c_abserr
        y[0] = r
        y[1] = z

        time = th
        new_time = th + time_step
        flag = 1
        flag = self.c_rkf45.r8_rkf45(y, yp, &time, new_time, &relerr, abserr, flag)
        return  y[0], y[1], time

    def run_n_steps(self, double r, double z, double th, double time_step, int n):
        """Perform n steps of RKF4(5) from point (r, z, th) in the Cylindrical
        coordinate system.

        Arguments:
            r (float): Radial position
            z (float): Vertical position
            th (float): Toroidal angle
            time_step (float): Size of time_step
            n (int): Number of steps of size time_step to perform

        Returns:
            r (float): New radial position
            z (float): Vertical position
            th (float): New toroidal angle
        """
        cdef:
            double y[2]
            double yp[2]
            double relerr, abserr, time, new_time
            int flag, i

        relerr = self.c_relerr
        abserr = self.c_abserr
        y[0] = r
        y[1] = z

        time = th
        new_time = th + time_step
        flag = 1
        for i in range(n):
            flag = self.c_rkf45.r8_rkf45(y, yp, &time, new_time, &relerr, abserr, flag)
            flag = 2
            new_time = time + time_step

        return  y[0], y[1], new_time