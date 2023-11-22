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

cdef class PyEmbreeAccell:

    cdef EmbreeAccell *c_eacc
    cdef dict name_to_mesh
    cdef list loaded_meshes_id
