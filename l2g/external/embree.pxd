# distutils: language = c++
# cython: language_level = 3

from libcpp cimport bool

cdef extern from "accell_embree.hpp" nogil:
    cdef cppclass EmbreeAccell:
        EmbreeAccell() except +
        void initializeDevice()
        unsigned int commitMesh(float* vertices, Py_ssize_t n_vertices,
                                unsigned* triangles, Py_ssize_t n_triangles)
        bool deleteMesh(unsigned geom_id)

        # RT functions
        void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                     double tnear, double tfar, int omp_thread)
        void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                     double tnear, double tfar)
        bool checkIfHit()
        bool checkIfHit(int omp_thread)

        int returnGeomId()
        int returnGeomId(int omp_thread)

        int returnPrimId()
        int returnPrimId(int omp_thread)

        void prepareThreadContainers() # Call this in serial before calling any
                                       # RT function
        void prepareThreadContainers(int num_threads)

cdef class PyEmbreeAccell:

    cdef EmbreeAccell *c_eacc
    cdef dict name_to_mesh
    cdef list loaded_meshes_id
