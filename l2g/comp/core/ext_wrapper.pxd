# distutils: language = c++
# cython: language_level = 3

# Wrapper class for both Cython/Python code

from libcpp cimport bool
from libcpp.vector cimport vector

from l2g.comp.core.ext_flt_cpp cimport FLT, EmbreeAccell

cdef class PyEmbreeAccell:

    cdef EmbreeAccell *c_eacc

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

    cdef  bool c_isHit(self) nogil
    cdef  bool c_isHit_omp(self, int omp_thread) nogil

    cdef  double c_getConlen(self) nogil
    cdef  double c_getConlen_omp(self, int omp_thread) nogil

    cdef  int c_getGeomID(self) nogil
    cdef  int c_getGeomID_omp(self, int omp_thread) nogil


    cdef  void getBCart(self, double r, double z, double phi, vector[double] &out) nogil
    cdef  void getBCyln(self, double r, double z, vector[double] &out) nogil
    cdef  double getPoloidalFlux(self, double r, double z) nogil
    cdef  double getFPol(self, double flux) nogil




