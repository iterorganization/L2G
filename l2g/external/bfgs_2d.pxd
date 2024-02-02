# distutils: language = c++
# cython: language_level = 3

from l2g.external.bicubic cimport BICUBIC_INTERP

cdef class PyBfgs2d:
    cdef BICUBIC_INTERP *c_bicubic
    cdef double bx1
    cdef double bx2
    cdef double by1
    cdef double by2
    cdef list _store

    cdef void c_set_interpolator(self, BICUBIC_INTERP *interp)

    cdef double line_search(self, double x, double y, double p1, double p2)

    # Zoom function for finding the adequate step length
    cdef double zoom(self, double x, double y, double a1, double a2, double f0, double gk0, double p1, double p2, double c1, double c2)
    # Interpolation minimizer functions
    cdef double quadmin(self, double a, double fa, double fpa, double b, double fb)
    cdef double cubicmin(self, double a, double fa, double fpa, double b, double fb, double c, double fc)