# distutils: language = c++
# cython: language_level = 3

from l2g.external.bicubic cimport BICUBIC_INTERP

cdef class PyBgfs2d:
    cdef BICUBIC_INTERP *c_bicubic

    cdef double line_search(self, double x, double y, double p1, double p2)
