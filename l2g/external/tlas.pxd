# distutils: language = c++
# cython: language_level = 3

from libcpp cimport bool

cdef extern from "tlas.hpp" nogil:
    cdef cppclass TLAS:
        TLAS() except +
        unsigned int commitMesh(float* vertices, Py_ssize_t n_vertices,
                                unsigned* triangles, Py_ssize_t n_triangles)
        bool deleteMesh(unsigned geom_id)

        # RT functions
        void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                     float tfar)

        bool checkIfHit()
        int returnGeomId()
        int returnPrimId()

cdef class PyTLAS:

    cdef TLAS *c_tlas
    cdef dict name_to_mesh
    cdef list loaded_meshes_id
