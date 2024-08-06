# distutils: language = c++
# cython: language_level = 3

from libcpp.vector cimport vector


cdef extern from "bicubic.hpp" nogil:
    struct BI_DATA:
        double r
        double z
        double val
        double valdx
        double valdy
        double dx
        double dy
        double minx
        double miny
        double maxx
        double maxy
        int nx
        int ny

    cdef cppclass BICUBIC_INTERP:


        BICUBIC_INTERP() except +

        void rcin(double x, double y, double &out_cell_x, double &out_cell_y)
        void interpolate(int r, int c)

        double m_dx, m_dy

        vector[double] mx, my
        int m_nx, m_ny

        vector[vector[double]] m_f, m_fdx, m_fdy, m_fdxdy

        void setArrays(vector[double] x, vector[double] y, vector[vector[double]] f)
        void populateContext(BI_DATA *context)
        void getValues(BI_DATA *context)
        void getSecondDerivativeValues(BI_DATA *context)

cdef class PyBicubic:
    cdef BICUBIC_INTERP *c_bicubic
    cdef BI_DATA c_BI_DATA
    cdef bint prepared
