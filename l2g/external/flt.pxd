# distutils: language = c++
# cython: language_level = 3

from l2g.external.embree cimport EmbreeAccell

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "flt.hpp" nogil:
    cdef enum FLT_OPTION "FLT::FLT_OPTION":
        BY_TIME=0
        BY_LENGTH=1

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
        void setIV(double r, double z, double phi)
        void setIV(double r, double z, double phi, int ompThread)
        void setTimeSpan(double tstop, double step)
        void setMaximumConnectionLength(double max_conlen)
        void setSelfIntersectionAvoidanceLength(double value)

        void setDirection(int direction)
        void setDirection(int direction, int ompThread)

        void setFltOption(int option)

        bool prepareInterpolation()
        void prepareThreadContainers() # Sequential. Prepares containers for 1
                                       # thread
        void prepareThreadContainers(int n) # For Openmp. Prepares containers
                                            # for n threads

        # Function for obtaining filedlines
        void getFL(vector[double] storage, bool with_flt) # Write to storage

        # Functions for running flt
        void runFLT() # Start FLT on given setIV
        void runFLT(int ompThread) # For running in OpenMP

        # Result containers
        vector[double] m_conlens
        vector[int] m_geom_hit_ids
        vector[int] m_prim_hit_ids

        void r8_flt(double t, double y[2])

        void getBCyln(double r, double z, vector[double] &out)
        void getBCart(double r, double z, double phi, vector[double] &out)
        double getPoloidalFlux(double r, double z)
        double getVacuumFPOL()
        void getPFValues(double r, double z, double &val, double &valdx, double &valdy, double &valdxdy, int ompThread)

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

    cdef bool flag_rarr, flag_zarr, flag_psi, flag_vfpol

    cpdef double getVacuumFPOL(self)
    cpdef void setShift(self, double rmove, double zmove)
    cpdef void setAbsError(self, double abserr)
    cpdef void setRelError(self, double relerr)

    cpdef void setIV(self, double r, double z, double phi)
    cdef  void c_setIV(self, double r, double z, double phi) nogil
    cdef  void c_setIV_omp(self, double r, double z, double phi, int omp_thread) nogil

    cpdef void setTimeSpan(self, double tstop, double step)
    cdef  void c_setTimeSpan(self, double tstop, double step) nogil

    cpdef void setMaximumConnectionLength(self, double max_conlen)
    cdef  void c_setMaximumConnectionLength(self, double max_conlen) nogil

    cpdef void setSelfIntersectionAvoidanceLength(self, double value)
    cdef  void c_setSelfIntersectionAvoidanceLength(self, double value) nogil

    cpdef void setDirection(self, int direction)
    cdef  void c_setDirection(self, int direction) nogil
    cdef  void c_setDirection_omp(self, int direction, int omp_thread) nogil

    cpdef void setFltOption(self, int option)

    cdef  void c_getFL(self, vector[double] &storage, bool with_flt) nogil

    cdef  void c_runFLT(self) nogil
    cdef  void c_runFLT_omp(self, int omp_thread) nogil

    cdef  void c_prepareThreadContainers(self, int omp_thread=*) nogil

    cdef  double c_getConlen(self) nogil
    cdef  double c_getConlen_omp(self, int omp_thread) nogil

    cdef  int c_getGeomID(self) nogil
    cdef  int c_getGeomID_omp(self, int omp_thread) nogil

    cdef int c_getPrimID(self) nogil
    cdef int c_getPrimID_omp(self, int omp_thread) nogil


    cdef  void getBCart(self, double r, double z, double phi, vector[double] &out) nogil
    cdef  void getBCyln(self, double r, double z, vector[double] &out) nogil
    cdef  double getPoloidalFlux(self, double r, double z) nogil