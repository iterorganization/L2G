# distutils: language = c++
# cython: language_level = 3

from libcpp cimport bool

cdef extern from "accell_embree.hpp" nogil:
    cdef cppclass EmbreeAccell:
        EmbreeAccell() except +
        unsigned int commitMesh(float* vertices, Py_ssize_t n_vertices,
                                unsigned* triangles, Py_ssize_t n_triangles)
        bool deleteMesh(unsigned geom_id)

        # RT functions
        void castRay(float ox, float oy, float oz, float dx, float dy, float dz,
                     double tnear, double tfar)

        bool checkIfHit()
        int returnGeomId()
        int returnPrimId()

cdef class PyEmbreeAccell:

    cdef EmbreeAccell *c_eacc
    cdef dict name_to_mesh
    cdef list loaded_meshes_id
