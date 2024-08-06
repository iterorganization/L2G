# distutils: language = c++
# cython: language_level = 3

from l2g.external.bicubic cimport BICUBIC_INTERP, BI_DATA

cdef class PyBfgs2d:
    cdef double bx1, bx2, by1, by2
    cdef list _store
    cdef BI_DATA c_BI_DATA
    # For some reason if the following line is above the doubles the code
    # randomly segfaults, when creating an instance of the class. Debug info
    # shows that it happens when assigning values to bx* and by*. It happens
    # when the code is compiled with -O3 and with gcc 11.3.1. Cannot reproduce
    # the problem with a minimal example...
    cdef BICUBIC_INTERP *c_bicubic

    cdef void c_set_interpolator(self, BICUBIC_INTERP *interp)

    cdef double line_search(self, double x, double y, double p1, double p2)

    # Zoom function for finding the adequate step length
    cdef double zoom(self, double x, double y, double a1, double a2, double f0, double gk0, double p1, double p2, double c1, double c2)
    # Interpolation minimizer functions
    cdef double quadmin(self, double a, double fa, double fpa, double b, double fb)
    cdef double cubicmin(self, double a, double fa, double fpa, double b, double fb, double c, double fc)
