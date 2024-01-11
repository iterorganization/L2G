# distutils: language = c++
# cython: language_level = 3

from l2g.external.bicubic cimport BICUBIC_INTERP

cdef class PyBgfs2d:
    cdef BICUBIC_INTERP *c_bicubic
    cdef double bx1
    cdef double bx2
    cdef double by1
    cdef double by2

    cdef double line_search(self, double x, double y, double p1, double p2)
    cdef void c_set_interpolator(self, BICUBIC_INTERP *interp)