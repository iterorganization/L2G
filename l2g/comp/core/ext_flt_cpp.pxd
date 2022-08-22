# distutils: language = c++
# cython: language_level = 3

from libcpp cimport bool
from libcpp.vector cimport vector

# External includes

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

        void r8_flt(double t, double y[2])

        void getBCyln(double r, double z, vector[double] &out)
        void getBCart(double r, double z, double phi, vector[double] &out)
        double getPoloidalFlux(double r, double z)
        double getVacuumFPOL()
        void debug_getValues(double r, double z, double &val, double &valdx, double &valdy, int ompThread)

        void setEmbreeObj(EmbreeAccell* accellObj)


cdef extern from "accell_embree.hpp" nogil:
    cdef cppclass EmbreeAccell:
        EmbreeAccell() except +
        void initializeDevice()
        unsigned int commitMesh(float* vertices, Py_ssize_t n_vertices,
                                unsigned* triangles, Py_ssize_t n_triangles)
        bool deleteMesh(unsigned geom_id)

