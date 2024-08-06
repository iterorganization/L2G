# distutils: language = c++
# cython: language_level = 3

from l2g.external.embree cimport EmbreeAccell

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "flt.hpp" nogil:
    cdef cppclass FLT:
        FLT() except +

        void setNDIM(Py_ssize_t NDIMR, Py_ssize_t NDIMZ)
        void setRARR(vector[double] r)
        void setZARR(vector[double] z)
        void setPSI(vector[double] psi)
        void setFARR(vector[double] psi) # FPOL domain
        void setFPOL(vector[double] fpol)
        void setVacuumFPOL(double vacuum_fpol)
        void setShift(double rmove, double zmove)

        void setAbsError(double abserr)
        void setRelError(double relerr)
        void setDesiredStep(double step)
        void setMaximumFieldlineLength(double max_fieldline_length)
        void setSelfIntersectionAvoidanceLength(double value)

        bool prepareInterpolation()
        void setNumberOfThreads(int n)
        void setPoints(vector[double] origin_points)
        void setStartingFLDirection(vector[int] directions)

        # Functions for running flt
        void runFLT() # Start FLT on given setIV

        # Function for obtaining filedlines
        void getFL(const double r, const double z, const double phi,
                   const int direction, vector[double] storage,
                   bool with_flt) # Write to storage

        # Result containers
        vector[double] m_out_fieldline_lengths
        vector[int] m_out_geom_hit_ids
        vector[int] m_out_prim_hit_ids

        void getBCyln(double r, double z, vector[double] &out)
        void getBCart(double r, double z, double phi, vector[double] &out)
        double getPoloidalFlux(double r, double z)
        double getVacuumFPOL()
        void setEmbreeObj(EmbreeAccell* accellObj)

# Function signatures that have c_ prefix are written so that they can be used
# in nogil blocks.

# All cdef functions are defined with nogil as the structure of the external
# program and calls are made to be thread safe, i.e., if there are threads
# involved, then the function signature passed the thread ID or functions
# that retrieve information from the interpolated data do not have any shared
# variables and the external libraries used are also thread safe.

# Function signatures that have the c_ prefix and _omp suffix are written so
# that they can be used in nogil parallel blocks, where threads have local
# containers.

cdef class PyFLT:

    cdef FLT *c_flt # C++ instance we are wrapping

