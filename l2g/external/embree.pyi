import numpy as np

class PyEmbreeAccell:
    """This class wraps the EmbreeAccell class from L2G_cpp and is used to
    provide geometries to Embree and holds information about the geometries,
    i.e., what was hit, etc...
    """
    def __init__(self) -> None: ...
    def commitMesh(self, vertices: np.ndarray, triangles: np.ndarray, name: str = "") -> int:
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
    def deleteMesh(self, geom_id: int) -> bool:
        """Removes a commited mesh from Embree if it exists

        Arguments:
            geom_id (unsigned): Non-negative ID of a geometry to remove.
        """
    def isEmpty(self) -> bool: ...
    def isMeshWithNameIn(self, name: str) -> bool: ...

    # RayTrace functions
    def castInfRay(self, ox: float, oy: float, oz :float, dx: float, dy: float, dz: float) -> None:
        """Cast an infinite ray from origin point (ox, oy, oz) with direction
        vector (dx, dy, dz).

        Arguments:
            ox (float): X-component of origin
            oy (float): Y-component of origin
            oz (float): Z-component of origin
            dx (float): X-component of direction
            dy (float): Y-component of direction
            dz (float): Z-component of direction
        """
    def castRay(self, ox: float, oy: float, oz :float, dx: float, dy: float, dz: float, tnear: float, tfar: float) -> None:
        """Cast a ray from origin point (ox, oy, oz) with direction
        vector (dx, dy, dz), starting at distance tnear and ending at distance
        tfar. tnear=0 means the ray starts at the (ox, oy, oz).

        Arguments:
            ox (float): X-component of origin
            oy (float): Y-component of origin
            oz (float): Z-component of origin
            dx (float): X-component of direction
            dy (float): Y-component of direction
            dz (float): Z-component of direction
            tfar (float): Distance along the direction.
        """
    # Result functions
    def checkIfHit(self) -> bool:
        """After calling castInfRay, check if there is a hit with this
        function.
        """
    # If checkIfHit == True, get cell ID number and loaded geometry number
    # with the following functions
    def returnGeomId(self) -> int:
        """Returns the id of the geometry that was hit. Geometry here
        corresponds to the mesh.
        """
    def returnPrimId(self) -> int:
        """Returns the id of the primitive that was hit. Here the id primitive
        corresponds to the cell (triangle) of the hit geometry.
        """
