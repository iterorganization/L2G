# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np

from libc.math cimport atan2
from libc.stdio cimport printf

from l2g.external.bicubic cimport PyBicubic
from l2g.external.rkf45 cimport PyRKF45FLT

cdef class PyRKF45FLT:
    def __cinit__(self):
        self.c_rkf45 = new RKF45()
        self.c_rkf45.set_omp_thread(0)

    def __init__(self):
        self.c_relerr = 1e-4
        self.c_abserr = 1e-4

    def setInterpolator(self, PyBicubic obj):
        self.c_rkf45.set_interpolator(obj.c_bicubic)

    def set_r_rmove(self, double r_move):
        self.c_rkf45.set_r_move(r_move)

    def set_z_rmove(self, double z_move):
        self.c_rkf45.set_z_move(z_move)

    def set_vacuum_fpol(self, double vacuum_fpol):
        self.c_rkf45.set_vacuum_fpol(vacuum_fpol)

    def run_step(self, double r, double z, double th, double time_step):

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