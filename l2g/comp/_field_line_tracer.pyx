# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np
import cython

## Cython imports
from l2g.external.flt cimport FLT
from l2g.external.embree cimport PyEmbreeAccell


# Std stuff
from libcpp cimport bool
from libc.stdio cimport printf
from libcpp.vector cimport vector
from libc.math cimport sqrt, atan2, acos

import os
import logging
from time import perf_counter
log = logging.getLogger(__name__)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int valSign(double inp) nogil:
    if inp > 0.0:
        return 1
    return -1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double angleBetweenVectors(vector[double] v1, vector[double] v2) noexcept nogil:
    cdef float dot, sq1, sq2, angle
    dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    sq1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]
    sq2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]
    angle = acos(dot / sqrt(sq1 * sq2))
    return angle

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void getBaryCenterInCyln(vector[double] points, float mul, vector[double] &out) noexcept nogil:
    """Get barycenter of triangle. And return it in cylindrical coordinate
    system.

    Arguments:
        points (np.float): 1D array of floats of size 9
        coord (str): Type of coordinates to return
        mul (float): Multiplier for dimension. For conversion from meter to mm,
                     or vice versa
        out (vector): Input/Output variable, here we store the (R, Z, Phi)
    """
    cdef double phi, r
    out[0] = (points[0] + points[3] + points[6]) / 3
    out[1] = (points[1] + points[4] + points[7]) / 3
    out[2] = (points[2] + points[5] + points[8]) / 3

    # Cylindrical
    phi = atan2(out[1], out[0])
    r = sqrt(out[0]**2 + out[1]**2)
    out[0] = mul * r
    out[1] = mul * out[2]
    out[2] = phi

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void getTriangleNormal(vector[double] points, vector[double] &out) nogil:
    """Get normal of triangle. Orientation of the points must be positive!

    Therefor the points list must be:
        [p1[0], p1[1], p1[2], p2[0]....]

    Arguments:
        points (np.float): 1D array of floats of size 9
    """
    cdef float x1,y1,z1,x2,y2,z2,norm

    x1 = points[3] - points[0]
    y1 = points[4] - points[1]
    z1 = points[5] - points[2]

    x2 = points[6] - points[0]
    y2 = points[7] - points[1]
    z2 = points[8] - points[2]

    out[0] = y1 * z2 - z1 * y2
    out[1] = z1 * x2 - x1 * z2
    out[2] = x1 * y2 - y1 * x2

    norm = sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2])
    out[0] = out[0] / norm
    out[1] = out[1] / norm
    out[2] = out[2] / norm


cdef class FieldLineTracer:
    """This is the class that performs FLT on input mesh data.

    Example usage:

    ``` py
    import l2g.comp

    flt = l2g.comp.FieldLineTracer()

    # Set the parameters, options and hlm settings

    flt.parameters
    flt.options
    flt.hlm_params

    # Set mesh data in the form of vertices and triangles
    flt.setTargetData(v, t)
    # Or just on points
    flt.setTargetData(v)

    # Set the equilibrium data
    flt.setEquilibrium(equilibrium)

    flt.applyParameters()
    flt.loadEq()
    flt.processMagneticData()

    # Run the FLT
    flt.runFLT()

    # Process the data
    flt.results
    flt.hlm_results

    # Get specific field lines
    flt.fl_ids = [1,2,3,...]

    flt.getFL()


    ```
    """
    cdef FLT *c_FLT
    cdef vector[double] c_points # In Cylindrical Coordinate system. (R [m], Z [m], Phi [rad])
    cdef vector[int] c_direction

    cdef public object name
    cdef public object hlm_params
    cdef public object parameters
    cdef public object options
    cdef public object embree_obj
    cdef public object equilibrium
    cdef public object eq
    # cdef public object point_results
    # cdef public object mesh_results
    cdef public object results
    cdef public object recalculate_magnetic_data
    cdef public object fl_results
    cdef public object owl_conlen_data
    cdef public object target_vertices
    cdef public object target_triangles
    cdef public object normals
    cdef public object fl_ids
    cdef public object hlm_results

    def __init__(self):
        self.name = "SFLT_case-1"
        import l2g.settings
        self.hlm_params: l2g.settings.HLM                 = l2g.settings.HLM()
        self.parameters: l2g.settings.Parameters          = l2g.settings.Parameters()
        self.options:    l2g.settings.Options             = l2g.settings.Options()
        import l2g.external.embree
        self.embree_obj: l2g.external.embree.PyEmbreeAccell = l2g.external.embree.PyEmbreeAccell()
        import l2g.equil
        self.equilibrium:l2g.equil.Equilibrium        = l2g.equil.Equilibrium()
        import l2g.external.equilibrium_analysis
        self.eq: l2g.external.equilibrium_analysis.EQA = l2g.external.equilibrium_analysis.EQA()

        #  Results objects
        # self.point_results: l2g.comp.L2GPointResults = l2g.comp.L2GPointResults() # Holds FLT result on points
        # self.mesh_results:  l2g.comp.L2GResults =      l2g.comp.L2GResults() # Holds FLT result on mesh
        import l2g.comp
        self.results: l2g.comp.L2GResults = l2g.comp.L2GResults()
        self.recalculate_magnetic_data: bool = True
        self.fl_results:    l2g.comp.L2GFLs =          l2g.comp.L2GFLs() # Holds fieldlines

        self.owl_conlen_data: np.ndarray = None

        # Geometry data for the target.
        self.target_vertices:  list = [] # Cartesian coordinate system. (X[m], Y[m], Z[m])
        self.target_triangles: list = []
        self.fl_ids:           list = []

        # HLM results
        self.hlm_results: l2g.comp.L2GResultsHLM = l2g.comp.L2GResultsHLM()

    def __cinit__(self):
        self.c_FLT = new FLT()

    def setParameters(self, parameters: 'l2g.comp.Parameters') -> None:
        self.parameters = parameters

    cpdef setEmbreeObj(self, PyEmbreeAccell embree_obj):
        self.embree_obj = embree_obj
        self.c_FLT.setEmbreeObj(embree_obj.c_eacc)

    def commitMeshesToEmbree(self, mesh_files: list[str] | str, dim_mul=1e-3) -> list:
        """Commits mesh from files to embree.

        Returns:
            ok (bool): Ok if finished

        Arguments:
            mesh_files: Either a single file or multiple files
            dim_mul (float): Dimension multiplier to convert mesh dimension to
                meters. Default 1e-3.
        """
        import l2g.mesh
        log.info("Commiting meshes to Embree object...")
        if isinstance(mesh_files, str):
            mesh_files = [mesh_files]

        out_ids = []

        for file in mesh_files:
            log.info(f"Commiting {file} to Embree object.")
            mesh = l2g.mesh.Mesh(file)
            v, t = mesh.getMeshData()
            out_ids.append(self.embree_obj.commitMesh(v * dim_mul, t))
            # Dimension is in meters

        return out_ids

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def setTargetData(self, vertices:np.ndarray,
                      triangles: np.ndarray | None = None) -> None:
        """Sets the target data.

        If vertices is provided and triangles is None, then the the vertices
        are already the origin points for the fieldlines.

        If vertices is provided and triangles is not None, then the vertices
        are the vertices of the triangles and triangles contains the ids of
        each triangle.

        The arrays can be whatever dimension as this function flattens it for
        further use.

        Additionally into the cpp vectors the actual origin points are stored,
        either the provided points or the barycenters. In cylindrical
        coordinate system!

        Arguments:
            vertices (np.ndarray): An array of points. These could already be
                the target points or are the vertices of the triangluar mesh.
                In meters
            triangles (np.ndarray | None): An optional numpy array of
                triangles. If this is not none, then the target data is
                considered a mesh. Otherwise the vertices are considered the
                target data.

        """

        if isinstance(vertices, list):
            vertices = np.asarray(vertices)
        if vertices.ndim > 1:
            vertices = vertices.flatten()
        if vertices.dtype != np.float32:
            vertices = vertices.astype(np.float32)
        self.target_vertices = vertices

        cdef:
            Py_ssize_t i3, i, n_tri
            double r, phi, dim_mul

            # Memory view on vertices
            float [:] m_vertices = vertices

        dim_mul = self.parameters.target_to_m

        self.c_points.clear()

        # Just points are provided. Copy them to the vector.
        if triangles is None:
            self.target_triangles = None
            self.normals = None
            # Target data
            # Copy the points to c_points

            ## Convert Cartesian points to Cylindrical points.
            r = self.target_vertices.shape[0] / 3.0
            # For some reason the following line raises an error in runtime
            # n_tri = <Py_ssize_t>(self.target_vertices.shape[0] / 3.0)
            n_tri = <Py_ssize_t> r
            self.c_points.resize(n_tri * 3)
            for i in range(n_tri):
                # Cylindrical
                i3 = 3*i

                phi = atan2(m_vertices[i3+1], m_vertices[i3])
                r = sqrt(m_vertices[i3]**2 + m_vertices[i3+1]**2)
                self.c_points[i3] = dim_mul * r # r
                self.c_points[i3+1] = dim_mul * m_vertices[i3+2] # z
                self.c_points[i3+2] = phi # phi
            # Only if there is a success with setting the point reset the
            # result object
            self.results.reset()
            return

        ## Mesh data! Some more steps to do
        if isinstance(triangles, list):
            triangles = np.asarray(triangles)
        if triangles.ndim > 1:
            triangles = triangles.flatten()

        self.target_triangles = triangles
        self.normals = np.empty(triangles.shape[0], np.float64)

        ## Get the Bary center
        cdef:
            Py_ssize_t pnt1, pnt2, pnt3
            vector[double] triangle_points
            vector[double] buff

            # Memory view on triangles and on vertices
            double [:] m_normals = self.normals

        self.c_points.resize(triangles.shape[0])
        buff.resize(3)
        triangle_points.resize(9)

        # TODO: try parallelizing this code! It is GIL friendly and should not
        # contain any overlapping things.
        for i in range(triangles.shape[0]/3):
            # Get the first point index
            i3 = 3 * i
            pnt1 = 3*triangles[i3]
            pnt2 = 3*triangles[i3+1]
            pnt3 = 3*triangles[i3+2]
            triangle_points[0] = m_vertices[pnt1]
            triangle_points[1] = m_vertices[pnt1+1]
            triangle_points[2] = m_vertices[pnt1+2]
            triangle_points[3] = m_vertices[pnt2]
            triangle_points[4] = m_vertices[pnt2+1]
            triangle_points[5] = m_vertices[pnt2+2]
            triangle_points[6] = m_vertices[pnt3]
            triangle_points[7] = m_vertices[pnt3+1]
            triangle_points[8] = m_vertices[pnt3+2]

            getBaryCenterInCyln(triangle_points, dim_mul, buff)
            # Copy it to the vector
            self.c_points[i3] = buff[0]
            self.c_points[i3+1] = buff[1]
            self.c_points[i3+2] = buff[2]
            # Copy it to the lsit attribute
            getTriangleNormal(triangle_points, buff)
            m_normals[i3] = buff[0]
            m_normals[i3+1] = buff[1]
            m_normals[i3+2] = buff[2]

        # Only if there is a success with setting the point reset the
        # result object
        self.results.reset()

    def setEquilibrium(self, equilibrium: 'l2g.equil.Equilibrium') -> None:
        """Set the equilibrium data, which is propagated to the eq analyze
        class and the external FLT C++ code.

        Arguments:
            equilibrium (l2g.equil.Equilibrium): Contains equilibrium magnetic
                data (2D axisymmetric)
        """
        self.equilibrium = equilibrium
        self.eq.setEquilibrium(equilibrium)
        self.parameters.plasma_r_displ = 0.0
        self.parameters.plasma_z_displ = 0.0
        self.recalculate_magnetic_data = True

    def applyParameters(self) -> None:
        """Propagates the parameters to the external FLT C++ code. Run this
        before running any FLT run functions.
        """
        log.info("Committing parameters")
        # Set the Embree pointer
        # self.c_FLT.applyRT(self.embree_obj)
        # Re-set the embree object. Should be checked.
        self.setEmbreeObj(self.embree_obj)

        # Apply RKF45 accuracy settings
        self.c_FLT.setAbsError(self.parameters.abs_error)
        self.c_FLT.setRelError(self.parameters.rel_error)

        # Apply plasma shift
        plasma_r_displ, plasma_z_displ = self.parameters.plasma_r_displ, self.parameters.plasma_z_displ
        self.c_FLT.setShift(plasma_r_displ, plasma_z_displ)
        # Apply plasma shift to the EQ object
        self.eq.setDisplacement(plasma_r_displ, plasma_z_displ)
        # Apply toroidal angle settings.
        self.c_FLT.setDesiredStep(self.parameters.time_step)
        # Apply maximum connection length settings
        self.c_FLT.setMaximumFieldlineLength(
            self.parameters.max_fieldline_length)
        # Apply the self intersection avoidance length.
        self.c_FLT.setSelfIntersectionAvoidanceLength(
            self.parameters.self_intersection_avoidance_length)

    def loadEq(self) -> None:
        """Loads the equilibrium data to the external FLT C++ code. Data is
        used for FL equations. Run this before running any FLT run functions.
        """
        log.info("Dumping equilibrium data to FLT object.")
        self.c_FLT.setNDIM(self.equilibrium.grid_dim_r,
                             self.equilibrium.grid_dim_z)
        self.c_FLT.setRARR(self.equilibrium.grid_r)
        self.c_FLT.setZARR(self.equilibrium.grid_z)
        self.c_FLT.setPSI(self.equilibrium.psi.flatten())
        self.c_FLT.setVacuumFPOL(self.equilibrium.fpol_vacuum)
        self.c_FLT.setFARR(self.equilibrium.fpol_flux)
        self.c_FLT.setFPOL(self.equilibrium.fpol)

        self.c_FLT.prepareInterpolation() # Prepares the interpolators

    def evaluateEq(self) -> None:
        """Tries to automatically determine the type of the equilibrium
        """
        log.info("Evaluating equilibrium")
        self.eq.evaluate()
        log.info(f"Type: {self.eq.getType()}, LCFS flux: {self.eq.getBoundaryFluxValue()} Webb/rad")

    @cython.boundscheck(False)
    def flipFieldLineDirections(self):
        cdef:
            int i

        for i in range(self.c_direction.size()):
            self.c_direction[i] = -self.c_direction[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def processMagneticData(self) -> None:
        """Here the magnetic data is processed. This function is called
        automatically before any fieldline tracing or fieldline obtaining
        function is called. But in any case if you need the data on the
        target points (either points in space or barycenter of triangular
        mesh) you can call this manually.
        """

        if not len(self.target_vertices):
            log.error("No target points loaded! Stopping")

        log.info("Processing magnetic data...")
        start = perf_counter()

        results = self.results
        cdef:
            Py_ssize_t n_points = self.c_points.size() // 3

        n_points = self.c_points.size() // 3
        if not self.recalculate_magnetic_data:
            log.info("Data was already calculated.")
            return
        self.recalculate_magnetic_data = False

        # # Allocate the arrays
        log.info("Allocating arrays")
        results.reset(n_points)
        log.info("Done")

        # First process the data that will be present in both cases.
        cdef:
            double [:] BVec = results.BVec
            double [:] BVecCyln = results.BVecCyln
            int [:] direction = results.direction
            double [:] flux = results.flux
            vector[double] buff1
            int vacuumFPolSign
            int i, i3

        if self.c_direction.size() != n_points:
            self.c_direction.resize(n_points)
        vacuumFPolSign = valSign(self.c_FLT.getVacuumFPOL())
        buff1.resize(3)
        log.info("Processing points")
        for i in range(n_points):
            i3 = 3*i
            self.c_FLT.getBCart(self.c_points[i3], self.c_points[i3+1],
                                self.c_points[i3+2], buff1)
            # Magnetic vector in Cartesian coordinate system
            BVec[i3] = buff1[0]
            BVec[i3+1] = buff1[1]
            BVec[i3+2] = buff1[2]
            self.c_FLT.getBCyln(self.c_points[i3], self.c_points[i3+1],
                                buff1)

            # Magnetic vector in Cylindrical coordinate system (tor, pol)
            BVecCyln[2*i] = buff1[0]
            BVecCyln[2*i+1] = buff1[1]

            # Magnetic poloidal flux [Web/rad]
            flux[i] = self.c_FLT.getPoloidalFlux(self.c_points[i3], self.c_points[i3+1])

            # Get the starting direction from vacuumFPolSign. Depending on
            # helicity or COCOS data notation, this is always -1 or positive
            # orientation
            direction[i] = vacuumFPolSign
            self.c_direction[i] = vacuumFPolSign

        ## Mesh data
        # Now check if we have mesh
        if self.target_triangles is None:
            # We do nothing
            results.empty = False
            log.info("Data processed")
            return

        # Now we continue. First memory views.
        cdef:
            # Normals are computed in the setTargetData!
            double [:] normals = self.normals
            double [:] angle = results.angle
            double [:] Bdot = results.Bdot
            vector[double] buff2
        buff2.resize(3)
        log.info("Obtaining direction for triangles.")
        for i in range(n_points):
            i3 = 3*i

            buff1[0] = BVec[i3]
            buff1[1] = BVec[i3+1]
            buff1[2] = BVec[i3+2]

            buff2[0] = normals[i3]
            buff2[1] = normals[i3+1]
            buff2[2] = normals[i3+2]

            # Get the angle
            angle[i] = angleBetweenVectors(buff1, buff2)

            # Get the Bdot
            Bdot[i] = buff1[0] * buff2[0] + buff1[1] * buff2[1] + buff1[2] * buff2[2]

            # Modify direction if necessary. If the normals is showing in the
            # opposite way of the magnetic vector or fieldline, we need to
            # follow the fieldline in the opposite direction.
            if Bdot[i] < 0:
                direction[i] = -direction[i]
                self.c_direction[i] = -self.c_direction[i]
        results.empty = False
        log.info("Data processed")

    def getFL(self) -> None:
        """Get fieldlines points. ONLY ON MESH DATA!
        """

        if not len(self.fl_ids):
            log.error("No target points to trace FL from. Stopping")
            return

        self.processMagneticData()

        # Just to be sure, convert the fl ids to numpy unsigned int 32
        self.fl_ids = np.asarray(self.fl_ids, np.uint32)

        # Make a check on the IDs, to see if the values are within interval
        # of triangles
        N_triangles = len(self.target_triangles)
        _stop = False
        for el in self.fl_ids:
            if el > N_triangles or el < 0:
                log.info(f"Skipping element id {el} as it points to no triangle")
                log.info("Please see the settings and remove invalid FL IDs.")
                _stop = True
        if _stop:
            return

        log.info("Getting FLs...")
        start = perf_counter()
        import l2g.comp

        self.fl_results = l2g.comp.L2GFLs()
        self.fl_results.target_to_m = self.parameters.target_to_m

        cdef:
            vector[vector[double]] c_fl_results
            vector[double] fl_buffer
            int i, tri_id, pnt_id

            # Memory views
            int [:] direction = self.results.direction

        for i in range(self.fl_ids.shape[0]):
            tri_id = self.fl_ids[i]
            pnt_id = 3 * tri_id
            fl_buffer.clear()
            self.c_FLT.getFL(self.c_points[pnt_id], self.c_points[pnt_id + 1],
                             self.c_points[pnt_id + 2],
                             direction[tri_id], fl_buffer,
                             self.options.switch_getFL_with_FLT)
            c_fl_results.push_back(fl_buffer)
        self.fl_results.points = <object> c_fl_results

        # Now scale the R, Z, depending on the target dim_mul, so that the
        # generated VTK has the same units as the target mesh.
        log.info(f"Finished getting FLs in {perf_counter() - start} seconds.")

    def getFLOnPoint(self, R: float, Z: float, Theta: float) -> list:
        """Obtain FL points that goes through the input parameters R, Z, Theta.
        """
        points = []

        cdef:
            vector[double] fl_buffer

        fl_buffer.clear()

        self.c_FLT.getFL(R, Z, Theta, 1, fl_buffer, self.options.switch_getFL_with_FLT)
        points.append(<object> fl_buffer)
        fl_buffer.clear()

        self.c_FLT.getFL(R, Z, Theta, -1, fl_buffer, self.options.switch_getFL_with_FLT)
        points.append(<object> fl_buffer)
        fl_buffer.clear()
        # Now let's order them.
        # points = points[0][::-1] + points[1]
        return points

    def runFLT(self) -> None:
        if not len(self.target_vertices):
            log.error("No target vertices loaded! Stopping")
            return

        if self.target_triangles is None:
            log.info("Starting FLT on points...")
        else:
            log.info("Starting FLT on mesh...")

        start = perf_counter()
        ## Setting the input data
        # Acquire the number of threads
        cdef int num_threads
        num_threads = os.cpu_count()
        if self.parameters.num_of_threads <= 0:
            user_num_threads = num_threads
        else:
            user_num_threads = min(self.parameters.num_of_threads, num_threads)
        self.c_FLT.setNumberOfThreads(user_num_threads)
        # Set the points from which FLs are tracked
        self.c_FLT.setPoints(self.c_points)
        # Set the starting direction of the points.
        self.c_FLT.setStartingFLDirection(self.c_direction)

        self.c_FLT.runFLT()

        # Copy the results!

        cdef int size, i
        size = self.c_FLT.m_out_fieldline_lengths.size()

        log.info("Copying results...")
        for i in range(size):
            self.results.conlen[i] = self.c_FLT.m_out_fieldline_lengths[i]
            self.results.geom_hit_ids[i] = self.c_FLT.m_out_geom_hit_ids[i]
            self.results.prim_hit_ids[i] = self.c_FLT.m_out_prim_hit_ids[i]
        log.info("Done")
        # Additionally create the mask which shows the final plasma pattern
        if not self.target_triangles is None:
            log.info("Creating wetted mask based on cutoff fieldline length")
            self.results.mask = self.applyShadowMask(
                np.ones(self.results.conlen.shape))

        log.info(f"Finished. It took {perf_counter() - start} seconds.")

    def obtainOuterWallConnectionGraph(self) -> None:
        """Obtains the outer wall connection length graph. Basically following
        FLs on the midplane.
        """

        radial_points = 500
        toroidal_points = 360

        points = self.createMidplanePoints(radialPoints=radial_points,
                                           toroidalSegments=toroidal_points)

        cdef:
            int i, size, fl_direction
            vector[double] c_points
            vector[int] c_direction
            double[:] m_conlen

        # First direction is down, so set the value accordingly!
        fl_direction = -1
        if self.c_FLT.getVacuumFPOL() < 0:
            fl_direction = 1

        size = len(points)
        c_points.resize(size)
        c_direction.resize(size, fl_direction)
        conlen = np.empty(size, np.float64)
        m_conlen = conlen # memoryview
        for i in range(size):
            c_points[i] = points[i]

        R_points = points[:3*radial_points:3]
        out = np.empty((radial_points, 4), np.float64)

        log.info("Obtaining OWL connection length graph")
        ## Getting the outer wall midplane profile.
        # First we obtain the outer wall - divertor fieldlines. This is what
        # is actually used by the ELM PLM.
        num_threads = os.cpu_count()
        if self.parameters.num_of_threads <= 0:
            user_num_threads = num_threads
        else:
            user_num_threads = min(self.parameters.num_of_threads, num_threads)
        self.c_FLT.setNumberOfThreads(user_num_threads)
        # Set the points from which FLs are tracked
        self.c_FLT.setPoints(c_points)
        # Set the starting direction of the points.
        self.c_FLT.setStartingFLDirection(c_direction)

        self.c_FLT.runFLT()

        # Copy the connection lengths values and flip the values of the
        # direction for connection lengths up.
        for i in range(size):
            # Connection length Down
            m_conlen[i] = self.c_FLT.m_out_fieldline_lengths[i]
            c_direction[i] = -c_direction[i]

        # Reshape them
        conlen_down = conlen.reshape((toroidal_points, radial_points))
        # Get the mean
        conlen_down_mean = np.mean(conlen_down, 0)
        # Store it
        out[:, 3] = conlen_down_mean

        # Fieldlines up!
        self.c_FLT.setStartingFLDirection(c_direction)
        self.c_FLT.runFLT()
        for i in range(size):
            # Connection length Up
            m_conlen[i] = self.c_FLT.m_out_fieldline_lengths[i]

        # Reshape them
        conlen_up = conlen.reshape((toroidal_points, radial_points))
        # Get the mean
        conlen_up_mean = np.mean(conlen_up, 0)
        # Store it
        out[:, 2] = conlen_up_mean

        psiLcfs = self.eq.getBoundaryFluxValue()
        Rb, _, _, Bpm = self.eq.getOWL_midplane()
        flux = self.point_results.flux[:radial_points]
        # Calculate drsep
        drsep = np.abs((flux - psiLcfs) / (Rb * Bpm))

        # Actual drsep. take from R_points
        drsep_r = R_points - R_points[0]

        out[:, 0] = drsep
        out[:, 1] = drsep_r
        self.owl_conlen_data = out

    def alignGeometryWithLCFS(self) -> None:
        """In most cases the wall silhouette used in generating the equilibrium
        do not coincide with the geometry used in FLT. Therefore in cases, such
        as limiter cases, where every millimeter counts, we require to align
        the input geometry, so that the equilibrium LCFS is limiting also on
        the geometry.

        !!! note

            In order to use it, the diagnostic for the equilibrium must have a
            the points of the LCFS or boundary contour prepared so that there
            are points to compare.

        The LCFS contour is clipped to the closest element on the mesh.
        Afterwards the distance is set as the displacement for the radial and
        vertical direction. If there are any original displacements, let's say
        vertical already applied, this is already taken into account since the
        LCFS points are following the shift, therefore the underlying C++ code
        will receive the full correct displacement.

        !!! todo

            Find a good algorithm to obtain LCFS contour points with a good
            resolution. But this has to be done in EQ class.

        """
        log.info("Checking whether the input geometry is aligned with the" +
                 " contact point of the data")

        # Simplest way to check is to check the Psi values. Note that this
        # will not work great, if for instance the poloidal flux map contains
        # a series of contour "islands" where the values overlap, meaning we
        # have a set of magnetic surfaces that have the same poloidal flux map.
        if self.results.empty:
            log.error("No FLT data.")
            return

        # Calculate the initial drsep.
        log.info("Initial check for drsep.")
        self.calculateDrsep()

        # Get the element that is closest.
        el_min = np.argmin(self.results.drsep)
        dr_min = self.results.drsep[el_min]
        elr_min = self.c_points[3 * el_min]
        elz_min = self.c_points[3 * el_min + 1]
        # psi_min = self.results.flux[el_min]
        log.info(f"Closest element {el_min}: R={elr_min} Z={elz_min}")
        log.info(f"Closest distance from boundary: {dr_min}")

        # Using the stored LCFS contour in the diagnostic class we find the
        # shortest path and use it as radial and vertical displacement.

        point = [elr_min, elz_min]

        # Calculate the distance**2 from the closest element to the closest
        # point on the LCFS contour.

        displ = self.eq.alignLcfsToPoint(point)
        if displ is None:
            log.info(f"Calculated displ is None. Stopping.")
            return
        log.info(f"Calculated displ = {displ}")

        prev_r_displ = self.parameters.plasma_r_displ
        prev_z_displ = self.parameters.plasma_z_displ
        self.parameters.plasma_r_displ = displ[0]
        self.parameters.plasma_z_displ = displ[1]
        log.info(f"Setting new plasma displacement to:")
        log.info(f"R={displ[0]}, Z={displ[1]}")
        # REAPPLY parameters
        self.applyParameters()
        # DO NOT CHANGE EQ data
        # With this we solidify that the LCFS is used from the input
        # equilibrium data and we shift the plasma only for the target/shadow
        # geometry. Otherwise if we move also the plasma in the input
        # equilibrium data then the LCFS changes again and we are in a cycle of
        # never aligning the mesh.
        self.eq.setDisplacement(prev_r_displ, prev_z_displ)
        self.loadEq()
        # Force recalculation of data!
        self.recalculate_magnetic_data = True
        self.processMagneticData()

    def calculateDrsep(self) -> None:
        """From the evaluated Flux data and the equilibrium data evaluate the
        radial distance along the midplane for the input target mesh data.
        """
        log.info("Calculating distance from the boundary on midplane for each mesh element.")

        # Evaluate the equilibrium
        self.eq.evaluate()
        # Get IWL or OWL parameters
        log.info(f"Calculating for side: {self.parameters.side}")

        Rb, _, _, Bpm = self.eq.getMidplaneInfo(which=self.parameters.side)
        boundary = self.eq.getBoundaryFluxValue()
        drsep = (self.results.flux - boundary) / (Rb * Bpm)

        # In case of some COCOS notations, the flux gradient direction can
        # go either away from the plasma or inside of the plasma. In other
        # words the the flux values are either multiplied with -1 or 1. Hence
        # the correction. Otherwise in case when the gradient goes inside the
        # drsep value will become negative as the flux values will be less then
        # the boundary value
        drsep *= self.equilibrium.psi_sign

        self.results.drsep = drsep

        eq_type = self.eq.getType()
        if eq_type == "div":
            # Also evaluate the distance from the 2nd separatrix.

            secondaryXPoint = self.eq.getSecondaryXFluxValue()

            if secondaryXPoint is not None:
                Rb, _, _, Bpm = self.eq.getMidplaneInfo(lcfs=secondaryXPoint)
                drsep2 = (self.results.flux - secondaryXPoint) / (Rb * Bpm)
                drsep2 *= self.equilibrium.psi_sign
                self.results.drsep2 = drsep2
            else:
                self.results.drsep2 = np.zeros(drsep.shape)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def createMidplanePoints(self, which:str='owl', length:float=0.3,
            _radialPoints:int=2000, toroidalSegments:int=360,
            angularStep:float=1) -> np.ndarray:
        """Based on the Equilibrium data, create points placed on the midplane
        plane, just outside the boundary contour (only first point should
        generally have connection lengths going to infinity, or free path on
        the boundary magnetic surface).

        Arguments:
            which (str): On which side, iwl = inner wall, owl = outer wall
                we wish to obtain the points
            length (float): Length in the major radius direction. That is,
                effectively to which drsep we wish to populate points. In
                meters.
            radialPoints (int): How many points in radial direction we wish to
                have.
            toroidalSegments (int): How many segments or lines of points do we
                wish to have in poloidal direction. The number of segments
                times the angularStep should at maximum be 360 degrees.
            angularStep (float): Angular distance between the segments. In
                degrees.

        Returns:
            points (array): 1D array of points: (R1, Z1, Phi1, R2, Z2, Phi2,...)
        """
        log.info("Creating points on the midplane")
        # First we need to obtain the R, Z point of the midplane, which sits
        # on the boundary magnetic surface.
        if self.eq is None:
            log.error("No EQ instance loaded! Stopping.")
            return

        self.eq.evaluate()

        # Obtain the R, Z point
        # Rb, Z, Btotal, Bpm = self.eq.get_midplane_info(which=which)
        _Rb, _Z, _, _ = self.eq.getMidplaneInfo()
        cdef:
            double Rb = _Rb
            double Z = _Z
            int radialPoints = _radialPoints
        log.info(f"Midplane height [m]: {Z}")
        log.info(f"Radial start [m]: {Rb}")
        log.info(f"Radial span [m]: {length}")
        log.info(f"Number of radial points per segment: {radialPoints}")
        log.info(f"Toroidal segments: {toroidalSegments}")
        log.info(f"Angular step: {angularStep}")
        log.info(f"Total angle: {angularStep * toroidalSegments}")

        cdef:
            int radial_dir, total_num_of_points, i, j, offset
            double radial_step, phi

        if which == 'owl':
            # Points go in +R direction away from owl midplane
            radial_dir = 1
        else:
            # Points go in -R direction away from iwl midplane
            radial_dir = -1

        # Now let's start creating the points.

        total_num_of_points = radialPoints * toroidalSegments * 3

        # Note: for post-processing, re-shaping will come in handy, so think of
        # segments as rows and radial points as columns
        # points.reshape((toroidalSegments, radialPoints * 3))

        points = np.empty(total_num_of_points, np.float32)
        cdef:
            float [:] m_points = points
        radial_step = length / (radialPoints - 1)
        phi = 0
        angularStep = np.deg2rad(angularStep) # Convert to radians
        # To ensure that the poloidal segments are rows
        for i in range(toroidalSegments):
            # Now to ensure radial points would be columns

            for j in range(radialPoints):
                r_point = Rb + radial_dir * j * radial_step

                # x = r_point * np.cos(phi)
                # y = r_point * np.sin(phi)
                offset = (i * radialPoints + j) * 3
                m_points[offset] = r_point
                m_points[offset + 1] = Z
                m_points[offset + 2] = phi

            phi += angularStep

        return points

    def applyHLM(self) -> None:
        """Applies an exponential plasma profile. Either single or double,
        depending on the type stored in self.hlm_params.
        """
        if self.equilibrium is None:
            log.error("No loaded equilibrium. Stopping")
            return

        if self.results.empty:
            log.error("No mesh results to apply Ramp Down on. Stopping")
            return

        if self.results.drsep is None:
            log.error("Mesh results are empty. Stopping.")
            return

        if not self.hlm_results.empty:
            self.hlm_results.reset()

        # Obtain the arrays from the FLT results
        drsep = self.results.drsep # In meters
        Bdot = np.abs(self.results.Bdot)
        n_cells = drsep.shape[0]
        bfield_mag = np.linalg.norm(self.results.BVec.reshape((n_cells, 3)), axis=1)
        self.eq.evaluate()

        Rb, _, Btotal, Bpm = self.eq.getMidplaneInfo(which=self.parameters.side)
        import l2g.hlm.general
        import l2g.hlm.steady_state
        import l2g.hlm.ramp_down

        if self.hlm_params.hlm_type == "single":
            q_par = l2g.hlm.general.single_exponential_psol(drsep=drsep,
                Bt=Btotal, Bpm=Bpm, Rb=Rb, P_sol=self.hlm_params.p_sol,
                F=0.5, lambda_q=self.hlm_params.lambda_q)
        elif self.hlm_params.hlm_type == "double":
            if self.hlm_params.longwave_misaligment_applied:
                q_par = l2g.hlm.general.double_exponential_psol_longwave(
                    drsep=drsep, Bt=Btotal, Bpm=Bpm, Rb=Rb,
                    lambda_q_main=self.hlm_params.lambda_q_main,
                    lambda_q_near=self.hlm_params.lambda_q_near,
                    Rq=self.hlm_params.ratio, P_sol=self.hlm_params.p_sol,
                    F=0.5, delta_mis=self.hlm_params.longwave_l)
            else:
                q_par = l2g.hlm.general.double_exponential_psol(drsep=drsep,
                    Bt=Btotal, Bpm=Bpm, Rb=Rb,
                    lambda_q_main=self.hlm_params.lambda_q_main,
                    lambda_q_near=self.hlm_params.lambda_q_near,
                    Rq=self.hlm_params.ratio, P_sol=self.hlm_params.p_sol, F=0.5)
        elif self.hlm_params.hlm_type == "custom":
            # We have points and profile
            q_par = l2g.hlm.general.custom(drsep=drsep,
                points=self.hlm_params.points, profile=self.hlm_params.profile,
                extrapolate=self.hlm_params.extrapolate,
                outside_value=self.hlm_params.outside_value)
            # The custom is defined on a midplane with magnetic surfaces going
            # perpendicular to it. In reality of course we have the toroidal
            # and poloidal component that already states that the flux tubes
            # go from the midplane at an angle. Therefore we multiply this
            # profile with the factor in order to actually scale it to the
            # midplane.
            # In other words apply the pitch at the midplane. Other functions
            # defined in the l2g.hlm already contains this term, except for the
            # q_parallel profile here.
            q_par *= Btotal / Bpm
        elif self.hlm_params.hlm_type == "elm":
            # The ELM data is loaded with an extra step.
            interELM = l2g.hlm.steady_state.inter_ELM(drsep, Rb, Btotal, Bpm,
                    Rb=self.hlm_params.r_break, P_sol=self.hlm_params.p_sol,
                    lambda_n=self.hlm_params.lambda_q_near,
                    lambda_m=self.hlm_params.lambda_q_main)
            elm = l2g.hlm.general.custom(drsep=drsep,
                points=self.hlm_params.points, profile=self.hlm_params.profile)
            q_par = interELM + elm

            Te = l2g.hlm.general.custom(drsep=drsep,
                points=self.hlm_params.points,
                profile=self.hlm_params.additional_profiles[0])

            Ti = l2g.hlm.general.custom(drsep=drsep,
                points=self.hlm_params.points,
                profile=self.hlm_params.additional_profiles[1])

            self.hlm_results.additional_arrays.append(elm)
            self.hlm_results.additional_arrays.append(interELM)
            self.hlm_results.additional_arrays.append(self.applyShadowMask(Te))
            self.hlm_results.additional_arrays.append(self.applyShadowMask(Ti))
        elif self.hlm_params.hlm_type == "l_mod":
            q_par = l2g.hlm.steady_state.inter_ELM(drsep=drsep,
                R_bdry=Rb, B_total=Btotal, B_pol=Bpm,
                Rb=self.hlm_params.r_break, P_sol=self.hlm_params.p_sol,
                lambda_n=self.hlm_params.lambda_q_near,
                lambda_m=self.hlm_params.lambda_q_main)

        elif self.hlm_params.hlm_type == "ramp-down":
            # P_sol taken from IMAS, hopefully.
            Ip = self.equilibrium.Ip
            if Ip < self.hlm_params.ip_transition:
                lambda_q = l2g.hlm.ramp_down.decay_length_L_mode_diverted(
                    a = self.equilibrium.a, Ip=Ip,
                    Area=self.equilibrium.Area, R=self.equilibrium.mag_axis_r)
                q_par = l2g.hlm.general.single_exponential_psol(drsep, Btotal,
                    Bpm, Rb, lambda_q * 1e-3, self.equilibrium.Psol)
            else:
                lambda_q = float("NaN")
                q_par = np.zeros(drsep.shape)

            self.hlm_results.additional_arrays = [lambda_q]
        else:
            q_par = np.zeros(drsep.shape)

        expansion = bfield_mag / Btotal
        q_inc = q_par * Bdot / Btotal

        q_inc = self.applyShadowMask(q_par * Bdot / Btotal)
        if self.hlm_params.hlm_type == "elm":
            # Also put incident arrays for the inter-ELM and ELM
            self.hlm_results.additional_arrays.append(self.applyShadowMask(elm * Bdot / Btotal))
            self.hlm_results.additional_arrays.append(self.applyShadowMask(interELM * Bdot / Btotal))

        self.hlm_results.q_inc = q_inc
        self.hlm_results.q_par = q_par
        self.hlm_results.flux_expansion = expansion

    def applyShadowMask(self, array: np.ndarray) -> np.ndarray:
        """This function applies the shadow mask to an input array.

        ONLY USED WITH MESH RESULTS!!!

        Several criteria are used for determining if a certain area is wetted:

         * Connection length of a field line:
           if conlen > cutoff = True else False
         * Geometries which mark fieldlines as shadowed. This is mainly used in
           regions where the FLs escape the chamber and are not captured by
           anything.

        Arguments:
            array (np.ndarray): A 1D array with N_cells elements on which the
                mask is applied

        Returns:
            masked_array (np.ndarray): A copy of the input array with the mask.

        """

        # First apply the connection length criteria

        out = None

        out = np.where(
            self.results.conlen >= self.parameters.cutoff_conlen,
            array, 0)

        # Now for every geom id in parameters.artificial_fl_catcher_geom_id
        # mark the fls as zero

        for geom_id in self.parameters.artificial_fl_catcher_geom_id:
            log.debug(f"Masking all elements with geom_hit_ids == {geom_id} to zero.")
            out = np.where(self.results.geom_hit_ids != geom_id,
                           out, 0)
        return out

    def debug_getCPoints(self):
        return <object> self.c_points
