# distutils: language = c++
# cython: language_level = 3

from libcpp.vector cimport vector
from l2g.external.bicubic cimport BICUBIC_INTERP

cdef extern from "rkf45.hpp" nogil:

    cdef cppclass RKF45:
        RKF45() except +

        void set_r_move(double r_move)
        void set_z_move(double z_move)
        void set_vacuum_fpol(double vacuum_fpol)

        void set_interpolator(BICUBIC_INTERP *interp)
        int r8_rkf45(double y[2], double yp[2], double *t, double tout,
                     double *relerr, double abserr, int flag)

cdef class PyRKF45FLT:
    cdef RKF45 *c_rkf45
    cdef double c_relerr
    cdef double c_abserr