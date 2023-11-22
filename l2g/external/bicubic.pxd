# distutils: language = c++
# cython: language_level = 3

from libcpp.vector cimport vector


cdef extern from "bicubic.hpp" nogil:

    cdef cppclass BICUBIC_INTERP:
        BICUBIC_INTERP() except +

        void rcin(double x, double y, double &out_cell_x, double &out_cell_y)
        void interpolate(int r, int c)

        double m_dx, m_dy

        vector[double] mx, my
        int m_nx, m_ny

        vector[vector[double]] m_f, m_fdx, m_fdy, m_fdxdy

        void prepareContainers(int number_of_omp_threads)
        void prepareContainers()

        void setArrays(vector[double] x, vector[double] y, vector[vector[double]] f)
        void getValues(double x, double y, double &val, double &valdx, double &valdy)
        void getAllValues(double x, double y, double &val, double &valdx, double &valdy, double &valdxdy)
        void getSecondDerivatives(double x, double y, double &valdxdx, double &valdydy)

cdef class PyBicubic:
    cdef BICUBIC_INTERP *c_bicubic
    cdef bint prepared
