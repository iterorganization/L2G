# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np

# External CPP class
from l2g.external.embree cimport EmbreeAccell

from libcpp.limits cimport numeric_limits

cdef class PyEmbreeAccell:
    """This class wraps the EmbreeAccell class from L2G_cpp and is used to
    provide geometries to Embree and holds information about the geometries,
    i.e., what was hit, etc...
    """

    def __cinit__(self):
        self.c_eacc = new EmbreeAccell()

    def __init__(self):
        # List of Loaded GeomIDs in the Embree object.
        self.loaded_meshes_id: list = []
        # Mapping of name to mesh if name is provided when supplying a mesh.
        # Useful when you want to dynamically see what meshes are loaded, to
        # load only when there is a new mesh. This is though entirely up to
        # the user to manage.
        self.name_to_mesh: dict = {}

    def commitMesh(self, vertices, triangles, name=""):
        """Commit a mesh, consisting of list of vertices and list
        of cells. The code reutrns the ID of the commited geometry, which the
        user then has to handle the managment.

        To avoid confusion when commiting meshes to EmbreeAccell the
        variables vertices and triangles can be normal python lists.
        Preferably the vertices should be a 1D array of floats (np.float32) and
        triangles should be a 1D array of unsigned integers (np.uint32).

        This functions handles:
            - 1D arrays that are not correct type. Remedied by casting.
            - 2D arrays that are of shape (3, N) (Column major) or (N, 3)
              (Row major), which is then flattened to (3*N) array.


        Arguments:
            vertices: Preferably 1D array of points. The size of the array
                      should be equal to 3 times the number of points.
            triangles: Preferably 1D array of indices of the triangles. The
                       size of the array should be equal to 3 times the
                       number of triangles.
            name (str, optional): Name of the geometry

        """
        cdef Py_ssize_t n_vertices, n_triangles
        cdef np.ndarray[float] npa_vertices
        cdef np.ndarray[unsigned int] npa_triangles

        # Process the input data in the correct form.
        try:
            if isinstance(vertices, list):
                vertices = np.array(vertices, np.float32)
            if vertices.dtype != np.float32:
                vertices = vertices.astype(np.float32)
            if len(vertices.shape) > 1:
                vertices = vertices.reshape(vertices.size)
        except Exception as e:
            print("Could not prepare vertices for committing to Embree.")
            raise e
        npa_vertices = vertices

        try:
            if isinstance(triangles, list):
                triangles = np.array(triangles, np.uint32)
            if triangles.dtype != np.uint32:
                triangles = triangles.astype(np.uint32)
            if len(triangles.shape) > 1:
                triangles = triangles.reshape(triangles.size)
        except Exception as e:
            print("Could not prepare triangles for committing to Embree.")
            raise e
        npa_triangles = triangles

        n_vertices = len(vertices)

        n_triangles = len(triangles) // 3

        geomId = self.c_eacc.commitMesh(&npa_vertices[0], n_vertices,
                                        &npa_triangles[0], n_triangles)
        self.loaded_meshes_id.append(geomId)
        if name != "":
            self.name_to_mesh[name] = geomId
        return geomId

    def deleteMesh(self, unsigned geom_id) -> bool:
        """Removes a commited mesh from Embree if it exists

        Arguments:
            geom_id (unsigned): Non-negative ID of a geometry to remove.
        """

        if not geom_id in self.loaded_meshes_id:
            return False
        ok = self.c_eacc.deleteMesh(geom_id)
        if ok:
            self.loaded_meshes_id.remove(geom_id)
            # Now to delete entry in the name_to_mesh.
            entry = [k for k, v in self.name_to_mesh.items() if v == geom_id]
            if entry:
                self.name_to_mesh.pop(entry[0])
        return ok

    def isEmpty(self) -> bool:
        if len(self.loaded_meshes_id):
            return False
        return True

    def isMeshWithNameIn(self, name: str) -> bool:
        if name in self.name_to_mesh:
            return True
        return False

    def castInfRay(self, float ox, float oy, float oz, float dx, float dy,
                   float dz):
        """Cast an infinite ray from origin point (ox, oy, oz) with direction
        vector (dx, dy, dz).
        """

        cdef:
            double tnear, tfar

        tnear = 0.05
        tfar = numeric_limits[double].infinity()
        self.c_eacc.castRay(ox, oy, oz, dx, dy, dz, tnear, tfar)

    def castRay(self, float ox, float oy, float oz, float dx, float dy,
                   float dz, float tnear, float tfar):
        """Cast an infinite ray from origin point (ox, oy, oz) with direction
        vector (dx, dy, dz).
        """

        self.c_eacc.castRay(ox, oy, oz, dx, dy, dz, tnear, tfar)

    def checkIfHit(self):
        """After calling castInfRay, check if there is a hit with this
        function.
        """
        return self.c_eacc.checkIfHit()

    def returnGeomId(self):
        """Returns the id of the geometry that was hit. Geometry here
        corresponds to the mesh.
        """
        return self.c_eacc.returnGeomId()

    def returnPrimId(self):
        """Returns the id of the primitive that was hit. Here the id primitive
        corresponds to the cell (triangle) of the hit geometry.
        """
        return self.c_eacc.returnPrimId()
