# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector

from l2g.comp.core.ext_flt_cpp cimport FLT, EmbreeAccell


cdef class PyEmbreeAccell:
    """This class wraps the EmbreeAccell class from L2G_cpp and is used to
    provide geometries to Embree and holds information about the geometries,
    i.e., what was hit, etc...



    """
    def __cinit__(self):
        self.c_eacc = new EmbreeAccell()

    def commitMesh(self, vertices, triangles, name=""):
        """Commit a mesh, consisting of list of vertices and list
        of cells. The code reutrns the ID of the commited geometry, which the
        user then has to handle the managment.

        To avoid confusion when commiting meshes to PyEmbreeAccell the
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
        return geomId

    def deleteMesh(self, unsigned geom_id):
        """Removes a commited mesh from Embree if it exists

        Arguments:
            geom_id (unsigned): Non-negative ID of a geometry to remove.
        """

        return self.c_eacc.deleteMesh(geom_id)

cdef class PyFLT:
    """This class wrapps the FLT class from L2G_cpp and houses functions, used
    either from Cython or Python. Hence some functions are "repeated" numerous
    times as they have a python declaration and Cython declaration.
    """

    def __cinit__(self):
        self.c_flt = new FLT()

    def setNDIM(self, Py_ssize_t NDIMR, Py_ssize_t NDIMZ):
        """Set the dimensions of the R,Z grid.

        Arguments:
            NDIMR (Py_ssize_t): Number of radial points
            NDIMZ (Py_ssize_t): Number of vertical points
        """
        # Py_ssize_t  corresponds to np.intp
        self.c_flt.setNDIM(NDIMR, NDIMZ)

    def setRARR(self, double[:] r):
        """Set the radial points array.

        Arguments:
            r (double[:]): 1D array, containing NDIMR points. In meters.
        """

        # Since the C++ class requires vectors, the data is copied. But since
        # the size is not so big, i.e., < 1000, this shouldn't cause
        # performance issues.
        cdef vector[double] _a
        for i in range(len(r)):
            _a.push_back(r[i])
        self.c_flt.setRARR(_a)

    def setZARR(self, double[:] z):
        """Set the vertical points array.

        Arguments:
            z (double[:]): 1D array, containing NDIMZ points. In meters.
        """
        # Since the C++ class requires vectors, the data is copied. But since
        # the size is not so big, i.e., < 1000, this shouldn't cause
        # performance issues.
        cdef vector[double] _a
        for i in range(len(z)):
            _a.push_back(z[i])
        self.c_flt.setZARR(_a)

    def setPSI(self, double[:] psi):
        """Sets the flux array. The array is provided as 1D, row oriented and
        the size must be equal to NDIMR * NDIMZ.

        Arguments:
            psi (double[:]): 1D array of poloidal magnetic flux values. In
                             Webb/rad.
        """

        # Since the C++ class requires vectors, the data is copied. But since
        # the size is not so big, i.e., < 10000, this shouldn't cause
        # performance issues.
        cdef vector[double] _a
        for i in range(len(psi)):
            _a.push_back(psi[i])
        self.c_flt.setPSI(_a)

    def setFARR(self, double[:] farr):
        """Sets the flux points, used for interpolating the FPol (poloidal
        current function) inside the plasma center. For FLT properties, only
        the FPol value in vacuum is required. So this function is not needed.

        Arguments:
            farr (double[:]): 1D array of flux points at which the FPol values
                              are defined. In Webb/rad.
        """

        # Since the C++ class requires vectors, the data is copied. But since
        # the size is not so big, i.e., < 1000, this shouldn't cause
        # performance issues.
        cdef vector[double] _a
        for i in range(len(farr)):
            _a.push_back(farr[i])
        self.c_flt.setFARR(_a)

    def setFPOL(self, double[:] fpol):
        """Sets the poloidal current function or FPol values defined over a
        series of flux points, set with above function. Again, same as with
        above function we do not require it, as this data is defined inside the
        plasma.

        Arguments:
            fpol (double[:]): 1D array of FPol values, defined over a series
                             of flux points. In m T.
        """

        # Since the C++ class requires vectors, the data is copied. But since
        # the size is not so big, i.e., < 1000, this shouldn't cause
        # performance issues.
        cdef vector[double] _a
        for i in range(len(fpol)):
            _a.push_back(fpol[i])
        self.c_flt.setFPOL(_a)

    def setVacuumFPOL(self, double fpol):
        """Sets the poloidal current function or FPol=Bt * R value in the
        vacuum. This value is used for calculating the toroidal component of
        the magnetic field.

        Arguments:
            fpol (double): Poloidal current function or FPol value in vacuum.
                           In m T (meter Tesla).
        """
        self.c_flt.setVacuumFPOL(fpol)

    cpdef double getVacuumFPOL(self):
        """Returns the FPol in Vacuum.
        Returns:
            vacuumFpol (double): FPol in Vacuum. In m T.
        """
        return self.c_flt.getVacuumFPOL()

    cpdef void setShift(self, double rmove, double zmove):
        """Sets the radial and vertical shift of the plasma.

        Arguments:
            rmove (double): Radial shift in meters.
            zmove (double): Vertical shift in meters.
        """
        self.c_flt.setShift(rmove, zmove)

    cpdef void setAbsError(self, double abserr):
        """Sets the absolute error. Less than 1. Default 1e-5

        Arguments:
            abserr (double): Absolute error used by RKF45.
        """
        self.c_flt.setAbsError(abserr)

    cpdef void setRelError(self, double relerr):
        """Sets the relative error. Less than 1. Default 1e-5

        Arguments:
            abserr (double): relative error used by RKF45.
        """
        self.c_flt.setRelError(relerr)

    cpdef void setIV(self, double r, double z, double phi):
        """Sets the initial point of a FL.

        Arguments:
            r (double): Radial position. In meters
            z (double): Vertical position. In meters
            phi (double): Toroidal angle position. In radians
        """
        self.c_flt.setIV(r, z, phi)

    cdef void c_setIV(self, double r, double z, double phi) nogil:
        """Same as setIV, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        self.c_flt.setIV(r, z, phi)

    cdef void c_setIV_omp(self, double r, double z, double phi, int omp_thread) nogil:
        """Same as setTimeSpan, except it is a
        Cython function, in order to use in Cython OpenMP parallel
        blocks.
        Arguments:
            omp_thread (int): Id of an OpenMP thread. Should be >= 0 or
                              < maximum number of OpenMP threads.
        """
        self.c_flt.setIV(r, z, phi, omp_thread)

    cpdef void setTimeSpan(self, double tstop, double step):
        """Sets the parametric time or toroidal angle end and the resolution
        at which we wish gather points or to check for intersections of a FL.

        Arguments:
            tstop (double): End of toroidal angle or parametric time. In
                            radians.
            step (double): Toroidal angle step at which we wish to gather
                           points. Used both for obtaining FL points and for
                           FLT intersection.
        """
        self.c_flt.setTimeSpan(tstop, step)

    cdef void c_setTimeSpan(self, double tstop, double step) nogil:
        """Same as setTimeSpan, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        self.c_flt.setTimeSpan(tstop, step)

    cpdef void setMaximumConnectionLength(self,double max_conlen):
        """Sets the maximum connection length to which we follow a FL.

        Arguments:
            max_conlen (double): Maximum connection length in meters.
        """
        self.c_flt.setMaximumConnectionLength(max_conlen)

    cdef void c_setMaximumConnectionLength(self,double max_conlen) nogil:
        """Same as setMaximumConnectionLength, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        self.c_flt.setMaximumConnectionLength(max_conlen)

    cpdef void setSelfIntersectionAvoidanceLength(self,double value):
        """Sets the initial length at which we avoid intersection tests. For
        example 0.005 (meters), this means that at the first 5 millimeters, do
        not check for intersection tests.

        Arguments:
            value (double): Initial length at which intersection tests are
                            ignored. In meters.
        """
        self.c_flt.setSelfIntersectionAvoidanceLength(value)

    cdef void c_setSelfIntersectionAvoidanceLength(self,double value) nogil:
        """Same as setSelfIntersectionAvoidanceLength, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        self.c_flt.setSelfIntersectionAvoidanceLength(value)

    cpdef void setDirection(self, int direction):
        """Sets the toroidal direction of FL tracing:
         - 1 - ACW (Anti ClockWise)
         - -1 - CW (ClockWise)

        The reason for this parameter is so that we always start tracing in the
        direction of the normals of the geometry. Hence this is the reason
        why the dot product between the magnetic field vector and the normal
        of a barycenter of a triangle is calculated and its sign checked.

        Arguments:
            direction (int): Toroidal direction of tracing.
        """
        self.c_flt.setDirection(direction)

    cdef void c_setDirection(self, int direction) nogil:
        """Same as setDirection, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        self.c_flt.setDirection(direction)

    cdef void c_setDirection_omp(self, int direction, int omp_thread) nogil:
        """Same as setDirection, except it is a
        Cython function, in order to use in Cython OpenMP parallel
        blocks.
        Arguments:
            omp_thread (int): Id of an OpenMP thread. Should be >= 0 or
                              < maximum number of OpenMP threads.
        """
        self.c_flt.setDirection(direction, omp_thread)

    cpdef void setFltOption(self, int option):
        """Sets the FLT option. Namely, the option tells the FLT which
        parameter is used to stop tracing a FL.

        Current options:
            0 - Stop tracing until the maximum set toroidal angle (tstop)
            1 - Stop tracing until the maximum connection length is achieved.

        Arguments:
            option (int):
        """
        self.c_flt.setFltOption(option)

    def prepare(self):
        """Calls the prepareInterpolation function of the FLT class.
        """
        cdef bool f
        f = self.c_flt.prepareInterpolation()
        return f

    def getFL(self, with_flt=False):
        """Python function for obtaining FL trajectory points.
        """
        cdef vector[double] storage
        self.c_flt.getFL(storage, with_flt)
        return <object> storage

    cdef void c_getFL(self, vector[double] &storage, bool with_flt) nogil:
        """Cython function (no Python GIL) for obtaining FL trajectory points.
        """
        self.c_flt.getFL(storage, with_flt)

    def runFLT(self):
        """Calls the runFLT function which starts the FLT.
        """
        self.c_flt.runFLT()

    cdef void c_runFLT(self) nogil:
        """Cython function (no Python GIL) for calling the runFLT function
        which starts the FLT.
        """
        self.c_flt.runFLT()

    cdef void c_runFLT_omp(self, int omp_thread) nogil:
        """Cython function (no Python GIL) for calling the runFLT function from
        an OpenMP parallel block.

        Arguments:
            omp_thread (int): Id of the OpenMP thread which calls the runFLT
                              function. The Id is used as an index for storing
                              thread local data, which is also managed by
                              the FLT C++ class.
        """
        self.c_flt.runFLT(omp_thread)

    def prepareThreadContainers(self):
        """Function that prepares the thread containers, albeit this is used in
        serial mode. This has to be called before calling runFLT.
        """
        self.c_flt.prepareThreadContainers()

    cdef void c_prepareThreadContainers(self, int omp_thread=1) nogil:
        """Cython function (no Python GIL) which calls the
        prepareThreadContainers of the C++ class. This allocates all the
        necessary vectors, storages and variables which is locally used by
        threads.

        No shared memory is used from the OpenMP side, instead this is done
        manually. This way when threads save its local data, it is saved in its
        allocated space, so that there is no thread clash.

        """
        self.c_flt.prepareThreadContainers(omp_thread)

    def isHit(self):
        """Returns the tracing status after following FL.

        This function is used directly from Python, or serial mode.

        Returns:
            hit (bool): If True, it means that the FL is intersected or
                        shadowed during tracing, else it isn't.
        """
        return self.c_flt.m_hits[0]

    cdef bool c_isHit(self) nogil:
        """Same as isHit, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        return self.c_flt.m_hits[0]

    cdef bool c_isHit_omp(self, int omp_thread) nogil:
        """Cython function (no Python GIL) for calling the runFLT function from
        an OpenMP parallel block.

        Returns:
            hit (bool): Returns the tracing result, performed by an OpenMP
                        thread with Id omp_thread.

        Arguments:
            omp_thread (int): Id of the OpenMP thread which calls the runFLT
                              function. The Id is used as an index for storing
                              thread local data, which is also managed by
                              the FLT C++ class.
        """
        return self.c_flt.m_hits[omp_thread]

    def getConLen(self):
        """Returns the connection length of a traced FL.

        This function is used directly from Python, or serial mode.

        Returns:
            conlen (double): Length of a FL.

        """
        return self.c_flt.m_conlens[0]

    cdef double c_getConlen(self) nogil:
        """Same as getConLen, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        return self.c_flt.m_conlens[0]

    cdef double c_getConlen_omp(self, int omp_thread) nogil:
        """Cython function (no Python GIL) for calling the runFLT function from
        an OpenMP parallel block.

        Returns:
            conlen (double): Length of a FL.

        Arguments:
            omp_thread (int): Id of the OpenMP thread which calls the runFLT
                              function. The Id is used as an index for storing
                              thread local data, which is also managed by
                              the FLT C++ class.
        """
        return self.c_flt.m_conlens[omp_thread]

    def getGeomID(self):
        """Gets the Id of the geometry with which the FL intersected.

        If this is called when no intersection happens, then you get the
        maximum value of an unsigned int.

        Returns:
            geomId (int): Id of intersected geometry.
        """
        return self.c_flt.m_geom_hit_ids[0]

    cdef int c_getGeomID(self) nogil:
        """Same as getGeomID, except it is a
        Cython function, in order to use in Cython code (no python GIL calls).
        """
        return self.c_flt.m_geom_hit_ids[0]

    cdef int c_getGeomID_omp(self, int omp_thread) nogil:
        """Cython function (no Python GIL) for calling the runFLT function from
        an OpenMP parallel block.

        Returns:
            geomId (int): Id of intersected geometry.

        Arguments:
            omp_thread (int): Id of the OpenMP thread which calls the runFLT
                              function. The Id is used as an index for storing
                              thread local data, which is also managed by
                              the FLT C++ class.
        """
        return self.c_flt.m_geom_hit_ids[omp_thread]

    cdef void getBCart(self, double r, double z, double phi, vector[double] &out) nogil:
        """Gets the Magnetic field vector in Cartesian coordinate system. Used
        when processing data on a geometry (getting quantities).

        Result is written in the out variable. (Bx, By, Bz)
        """
        self.c_flt.getBCart(r, z, phi, out)

    cdef void getBCyln(self, double r, double z, vector[double] &out) nogil:
        """Gets the Magnetic field vector (poloidal and toroidal component)
        in Cylindrical coordinate system. Used when processing data on a
        geometry (getting quantities).

        Result is written in the out variable. (Bpol, Btor)

        Arguments:
            r (double): Radial position
            z (double): Vertical position
            out (vector, inout): Vector for writing values
        """
        self.c_flt.getBCyln(r, z, out)

    cdef double getPoloidalFlux(self, double r, double z) nogil:
        """Gets the polodidal magnetic flux value. Used when processing data on
        a geometry (getting quantities).

        Returns:
            flux (double): Poloidal magnetic flux. In Webb/rad.

        Arguments:
            r (double): Radial position
            z (double): Vertical position
        """
        cdef double value
        value = self.c_flt.getPoloidalFlux(r, z)
        return value

    cdef double getFPol(self, double flux) nogil:
        """Gets the poloidal current function (FPol = Bt R) value. Used when
        processing data on a geometry (getting quantities).

        Returns:
            fpol (double): FPol value. In m T.

        Arguments:
            flux (double): Magnetic poloidal flux value at which we wish to
                           obtain FPol.
        """
        cdef double value
        value = self.c_flt.getFPol(flux)
        return value

    def applyRT(self, PyEmbreeAccell obj):
        """This function is used to bind an EmbreeAccell C++ object to the
        FLT C++ object. Only pointer is used, so essentially we could switch
        different Embree objects.
        """
        self.c_flt.setEmbreeObj(obj.c_eacc)
