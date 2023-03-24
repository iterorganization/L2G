# distutils: language = c++
# cython: language_level = 3

# This file contains the core function that are exposed to python for running
# FLT functions.

# Import the wrapped external classed
from l2g.comp.core.ext_wrapper cimport PyEmbreeAccell, PyFLT
from l2g.comp import L2GResults, L2GPointResults, L2GFLs


import numpy as np
cimport numpy as np

# Std stuff
from libcpp cimport bool
from libc.stdio cimport printf
from libcpp.vector cimport vector
from libc.math cimport sqrt, atan2, acos

from cython.parallel import prange, threadid, parallel
import cython


cimport openmp # For getting number of threads

cdef int get_num_threads():
    cdef int num_threads
    with nogil, parallel():
        num_threads = openmp.omp_get_num_threads()
        with gil:
            return num_threads

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void getBaryCenter(vector[double] points, float mul, vector[double] &out) nogil:
    """Get barycenter of triangle.

    Arguments:
        points (np.float): 1D array of floats of size 9
        coord (str): Type of coordinates to return
        mul (float): Multiplier for dimension. For conversion from meter to mm,
                     or vice versa
        out (vector): Input/Output variable, here we store the (R, Z, Phi)
    """
    cdef float phi, r
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double angleBetweenVectors(vector[double] v1, vector[double] v2) nogil:
    cdef float dot, sq1, sq2, angle
    dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    sq1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]
    sq2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]
    angle = acos(dot / sqrt(sq1 * sq2))
    return angle

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int valSign(double inp) nogil:
    if inp > 0.0:
        return 1
    return -1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef processDataOnPoints(PyFLT flt_obj, float [:] tg_vertices,
                          results, double dim_mul=1e-3):
    """Function that pre-precoess input target points data. From the data and
    the loaded equilibrium, quantities relevant to FLT analysis are calculated:

    * BVec
    * flux
    * direction
    """
    cdef:
        Py_ssize_t N_vertices = tg_vertices.shape[0] // 3
        Py_ssize_t i, p1, p2, p3, iS
        double t, r, phi

    # Variable results should be provided through python.
    # results = L2GPointResults()

    results.points = np.empty(N_vertices * 3, np.float)
    results.BVec = np.empty(N_vertices * 3, np.float)
    results.direction = np.empty(N_vertices, np.int32)
    results.flux = np.empty(N_vertices, np.float)

    # If results haven't been made, then populate mask and conlen with
    # empty arrays or zeroed arays
    results.conlenUp = np.empty(N_vertices, dtype=np.float)
    results.conlenDown = np.empty(N_vertices, dtype=np.float)

    results.geom_hit_ids_up = np.empty(N_vertices, dtype=np.int32)
    results.prim_hit_ids_up = np.empty(N_vertices, dtype=np.int32)
    results.geom_hit_ids_down = np.empty(N_vertices, dtype=np.int32)
    results.prim_hit_ids_down = np.empty(N_vertices, dtype=np.int32)

    cdef:
        double [:] points = results.points
        double [:] BVec = results.BVec
        double [:] flux = results.flux
        int [:] direction = results.direction

    cdef:
        vector[double] vecBuff
        vector[double] vecBuff2
        int vacuumFPolSign
    vacuumFPolSign = valSign(flt_obj.getVacuumFPOL())
    vecBuff.resize(3)
    vecBuff2.resize(2)
    # For naive implementation this is required!
    flt_obj.c_prepareThreadContainers()
    for i in range(N_vertices):
        iS = 3 * i

        flt_obj.getBCart(tg_vertices[iS], tg_vertices[iS + 1],
                         tg_vertices[iS + 2], vecBuff)
        BVec[iS] = vecBuff[0]
        BVec[iS + 1] = vecBuff[1]
        BVec[iS + 2] = vecBuff[2]
        # For now no use for the following?
        # flt_obj.getBCyln(tg_vertices[iS], tg_vertices[iS + 1], vecBuff2)
        # BVecCyln[2*i] = vecBuff2[0]
        # BVecCyln[2*i + 1] = vecBuff2[1]

        # Determine the direction for Up and Downward. Same procedure as in FLT
        # except here we have standard normals (0, 0, -1) and (0, 0, 1)
        # So, we need only to check the Z component of the magnetic vector and
        # that is it.
        direction[i] = vacuumFPolSign
        if vecBuff[2] < 0.0:
            direction[i] = -vacuumFPolSign

        # # Transform point into R, Z, Phi. Well just R and Phi
        # r = sqrt(tg_vertices[iS]**2 + tg_vertices[iS + 1]**2) * dim_mul
        # phi = atan2(tg_vertices[iS + 1], tg_vertices[iS])

        # Store the values inside points
        points[iS] = tg_vertices[iS] # R
        points[iS + 1] = tg_vertices[iS + 1] # Z
        points[iS + 2] = tg_vertices[iS + 2] # Phi
        flux[i] = flt_obj.getPoloidalFlux(tg_vertices[iS], tg_vertices[iS + 1])
    results.empty = False

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef processData(PyFLT flt_obj, float [:] tg_vertices,
                  unsigned int [:] tg_cells, results,
                  double dim_mul=1e-3):
    """Function that pre-process input target triangular mesh data. From the
    data and the loaded equilibrium, quantities relevant to FLT analysis are
    calculated:

     * Triangle barycenters in (R, Z, Phi) space of units (m, m, rad)
     * Triangle normals
     * BVec: Magnetic field vector in Cartesian coordinate system on barycenters in T
     * BDot: Dot product between magnetic field vector and normals on barycenters in T
     * flux: Poloidal flux value in Webb/rad on barycenters
     * angle: Angle between the magnetic field vector and normal on barycenter
     * direction: The direction of the FL on barycenter

    Arguments:
        flt_obj (PyFLT): FLT object that holds information about equilibrium
        tg_vertices (list): 1D array of floats containing target points
        tg_cells (list): 1D array of uint32 containing triangles.
        dim_mul (double): Dimension multiplier that transform input target data
            to meters.
    Returns:
        results (L2GResult): Result class populated with data.
    """

    cdef:
        Py_ssize_t N_vertices = tg_vertices.shape[0] // 3
        Py_ssize_t N_cells = tg_cells.shape[0] // 3
        Py_ssize_t i, p1, p2, p3, iS
        double t

    # Variable results must be provided through python. Cython code should only
    # be used to quicken parts of the code that would benefit from it.
    # results = L2GResults()
    results.baryCent = np.empty(N_cells * 3, np.float)
    results.normals = np.empty(N_cells * 3, np.float)
    results.BVec = np.empty(N_cells * 3, np.float)
    results.BVecCyln = np.empty(N_cells * 2, np.float)
    results.Bdot = np.empty(N_cells, np.float)
    results.flux = np.empty(N_cells, np.float)
    results.angle = np.empty(N_cells, np.float)
    results.direction = np.empty(N_cells, np.int32)
    results.geom_hit_ids = np.empty(N_cells, np.int32)
    results.prim_hit_ids = np.empty(N_cells, np.int32)

    # Assign the target_verices and target_cells to results
    results.vertices = tg_vertices
    results.triangles = tg_cells

    # If results haven't been made, then populate mask and conlen with
    # empty arrays or zeroed arays
    results.mask = np.empty(N_cells, dtype=np.float)
    results.conlen = np.empty(N_cells, dtype=np.float)

    cdef:
        double [:] baryCent = results.baryCent
        double [:] normals = results.normals
        double [:] BVec = results.BVec
        double [:] BVecCyln = results.BVecCyln
        double [:] Bdot = results.Bdot
        double [:] flux = results.flux
        double [:] angle = results.angle
        int [:] direction = results.direction

    # Required buffers
    cdef:
        vector[double] vecBuff1 # For B_vec in Cartesian coordinate system
        vector[double] vecBuff2 # For triangle normal
        vector[double] vecBuff3 # For (B_poloidal, B_toroidal)
        vector[double] trianglePoints
        int vacuumFPolSign
    trianglePoints.resize(9)
    vecBuff1.resize(3)
    vecBuff2.resize(3)
    vecBuff3.resize(2)
    vacuumFPolSign = valSign(flt_obj.getVacuumFPOL())
    # For naive implementation this is required!
    flt_obj.c_prepareThreadContainers()
    for i in range(N_cells):
        iS = 3 * i
        p1 = 3 * tg_cells[i * 3]
        p2 = 3 * tg_cells[i * 3 + 1]
        p3 = 3 * tg_cells[i * 3 + 2]


        trianglePoints[0] = tg_vertices[p1]
        trianglePoints[1] = tg_vertices[p1 + 1]
        trianglePoints[2] = tg_vertices[p1 + 2]

        trianglePoints[3] = tg_vertices[p2]
        trianglePoints[4] = tg_vertices[p2 + 1]
        trianglePoints[5] = tg_vertices[p2 + 2]

        trianglePoints[6] = tg_vertices[p3]
        trianglePoints[7] = tg_vertices[p3 + 1]
        trianglePoints[8] = tg_vertices[p3 + 2]
        # vecBuff = getBaryCenter(trianglePoints, mul=1e-3)
        getBaryCenter(trianglePoints, dim_mul, vecBuff1)
        baryCent[iS] = vecBuff1[0]
        baryCent[iS + 1] = vecBuff1[1]
        baryCent[iS + 2] = vecBuff1[2]
        flt_obj.getBCart(baryCent[iS], baryCent[iS + 1], baryCent[iS + 2], vecBuff1)
        BVec[iS] = vecBuff1[0]
        BVec[iS + 1] = vecBuff1[1]
        BVec[iS + 2] = vecBuff1[2]
        flt_obj.getBCyln(baryCent[iS], baryCent[iS + 1], vecBuff3)
        BVecCyln[2*i] = vecBuff3[0]
        BVecCyln[2*i + 1] = vecBuff3[1]
        getTriangleNormal(trianglePoints, vecBuff2)
        normals[iS] = vecBuff2[0]
        normals[iS + 1] = vecBuff2[1]
        normals[iS + 2] = vecBuff2[2]
        flux[i] = flt_obj.getPoloidalFlux(baryCent[iS], baryCent[iS + 1])
        Bdot[i] = vecBuff1[0] * vecBuff2[0] + vecBuff1[1] * vecBuff2[1] + vecBuff1[2] * vecBuff2[2]
        angle[i] = angleBetweenVectors(vecBuff1, vecBuff2)

        direction[i] = vacuumFPolSign
        if Bdot[i] < 0:
            direction[i] = -direction[i]
    results.empty = False

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef runFLT(PyFLT flt_obj, float [:] tg_vertices,
             unsigned int [:] tg_cells, int user_num_threads=1,
             results=None, double dim_mul=1e-3):
    """Runs the FLT trace on the provided target data.

    All results are saved into the results python object. It is a simple class
    containing variables, where to store data.

    Arguments:
        flt_obj (PyFLT): The FLT objects
        tg_vertices (list): 1D array of points.
        tg_cells (list): 1D array of triangles.
        user_num_threads (int): Number of OpenMP threads to use
        results (L2GResults): container for results. If empty a new one is made
        dim_mul (float): Dimension multipler that transfroms target to meters!

    Returns:
        results (L2GResults): Class that holds FLT results

    """

    cdef:
        Py_ssize_t N_vertices = tg_vertices.shape[0] // 3
        Py_ssize_t N_cells = tg_cells.shape[0] // 3
        Py_ssize_t i, iS
        Py_ssize_t threadId
        int num_threads
        double t

    # Acquire the number of threads
    num_threads = get_num_threads()

    if user_num_threads == 0:
        user_num_threads = num_threads
    elif user_num_threads > num_threads:
        user_num_threads = num_threads
    # Prepare Results array

    # Pre-processing, i.e., obtaining all necessary values that are used in the
    # FLT or applying HLM must be done outside this function.
    # if results is None:
        # results = processData(flt_obj, tg_vertices, tg_cells, dim_mul)

    # Memory view on data
    cdef:
        double [:] baryCent = results.baryCent
        double [:] normals = results.normals
        double [:] BVec = results.BVec
        double [:] Bdot = results.Bdot
        double [:] flux = results.flux
        double [:] angle = results.angle
        double [:] mask = results.mask
        double [:] conlen = results.conlen
        int [:] direction = results.direction
        int [:] geom_hit_ids = results.geom_hit_ids
        int [:] prim_hit_ids = results.prim_hit_ids
        double [:] initial_value

    flt_obj.c_prepareThreadContainers(user_num_threads)
    # Parallel work

    # Use barycenters of triangles as starting point
    initial_value = baryCent

    with nogil, parallel(num_threads=user_num_threads):
        threadId = threadid()
        for i in prange(N_cells, schedule="dynamic"):
            # threadId = threadid()
            iS = 3 * i
            flt_obj.c_setDirection_omp(direction[i], threadId)
            flt_obj.c_setIV_omp(initial_value[iS], initial_value[iS + 1], initial_value[iS + 2],
                               threadId)

            flt_obj.c_runFLT_omp(threadId)
            conlen[i] = flt_obj.c_getConlen_omp(threadId)
            geom_hit_ids[i] = flt_obj.c_getGeomID_omp(threadId)
            prim_hit_ids[i] = flt_obj.c_getPrimID_omp(threadId)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef runFLTonPoints(PyFLT flt_obj, float [:] tg_vertices, int user_num_threads=1,
                     results=None, double dim_mul=1e-3):
    """Runs the FLT trace on the provided target data (points).

    Essentially this could be used for obtaining connection length graphs on
    midplane, since we do not propagate information about normals or angles.

    TODO, maybe: Include normals, as (0, 0, 1) and (0, 0, -1). If Midplanes are
    targeted then angles should be perpendicular. To implement in
    prepareDataOnPoints

    Arguments:
        flt_obj (PyFLT): The FLT object
        tg_vertices (list): 1D array of points.
        user_num_threads (int): Number of OpenMP threads to use
        results (L2GPointResults): container for results. If empty a new one is made
        dim_mul (float): Dimension multipler that transfroms target to meters!

    Returns:
        results (L2GPointResults): Class that holds FLT results
    """
    cdef:
        Py_ssize_t N_vertices = tg_vertices.shape[0] // 3
        Py_ssize_t i, iS, direction
        Py_ssize_t threadId
        int num_threads
        double t

    # Acquire the number of threads
    num_threads = get_num_threads()

    if user_num_threads == 0:
        user_num_threads = num_threads
    elif user_num_threads > num_threads:
        user_num_threads = num_threads
    # Prepare Results array

    if results is None:
        results = processDataOnPoints(flt_obj, tg_vertices, dim_mul)

    # Memory view on data
    cdef:
        double [:] points = results.points
        double [:] BVec = results.BVec
        double [:] flux = results.flux
        double [:] conlenUp = results.conlenUp
        double [:] conlenDown = results.conlenDown
        int [:] direction_map = results.direction
        int [:] geom_hit_ids_up = results.geom_hit_ids_up
        int [:] prim_hit_ids_up = results.prim_hit_ids_up
        int [:] geom_hit_ids_down = results.geom_hit_ids_down
        int [:] prim_hit_ids_down = results.prim_hit_ids_down

    flt_obj.c_prepareThreadContainers(user_num_threads)

    with nogil, parallel(num_threads=user_num_threads):
        threadId = threadid()

        for i in prange(N_vertices, schedule="dynamic"):
            # First go Down!
            direction = -1 * direction_map[i]
            iS = 3 * i
            flt_obj.c_setDirection_omp(direction, threadId)
            flt_obj.c_setIV_omp(points[iS], points[iS + 1], points[iS + 2],
                               threadId)
            flt_obj.c_runFLT_omp(threadId)
            conlenDown[i] = flt_obj.c_getConlen_omp(threadId)
            geom_hit_ids_down[i] = flt_obj.c_getGeomID_omp(threadId)
            prim_hit_ids_down[i] = flt_obj.c_getPrimID_omp(threadId)
            # Next go Up!
            direction = direction_map[i]
            flt_obj.c_setDirection_omp(direction, threadId)
            flt_obj.c_setIV_omp(points[iS], points[iS + 1], points[iS + 2],
                               threadId)

            flt_obj.c_runFLT_omp(threadId)
            conlenUp[i] = flt_obj.c_getConlen_omp(threadId)
            geom_hit_ids_up[i] = flt_obj.c_getGeomID_omp(threadId)
            prim_hit_ids_up[i] = flt_obj.c_getPrimID_omp(threadId)
    return results

def getFL(PyFLT flt_obj, float [:] tg_vertices, unsigned int [:] tg_cells,
          unsigned int [:] FL_ids, results=None, double dim_mul=1e-3,
          bool with_flt=False):
    """Gets FLs based on selected element IDs.

    Arguments:
        flt_obj (PyFLT): The FLT objects
        tg_vertices (list): 1D array of points.
        FL_ids (list): 1D array of triangles.
        user_num_threads (int): Number of OpenMP threads to use
        results (L2GResults): container for results. If empty a new one is made
        dim_mul (float): Dimension multipler that transfroms target to meters!

    Returns:
        results (L2GResults): Class that holds FLT results

    """

    cdef:
        Py_ssize_t N_vertices = tg_vertices.shape[0] // 3
        Py_ssize_t N_cells = FL_ids.shape[0]
        Py_ssize_t i, iS
        double t

    # Prepare Results array

    if results is None:
        results = processData(flt_obj, tg_vertices, tg_cells, dim_mul)

    fl_results = L2GFLs()

    cdef vector[vector[double]] c_fl_results

    # Memory view on data
    cdef:
        double [:] baryCent = results.baryCent
        double [:] normals = results.normals
        double [:] BVec = results.BVec
        double [:] Bdot = results.Bdot
        double [:] flux = results.flux
        double [:] angle = results.angle
        double [:] mask = results.mask
        double [:] conlen = results.conlen
        int [:] direction = results.direction

    cdef unsigned int el
    cdef vector[double] fl_buffer
    # Prepares containers for one thread
    flt_obj.c_prepareThreadContainers()
    # Trace limited by time

    for i in range(N_cells):
        el = FL_ids[i]
        iS = 3 * el
        flt_obj.c_setDirection(direction[el])
        flt_obj.c_setIV(baryCent[iS], baryCent[iS + 1], baryCent[iS + 2])
        fl_buffer.clear()
        flt_obj.c_getFL(fl_buffer, with_flt)

        c_fl_results.push_back(fl_buffer)
    fl_results.points = <object> c_fl_results
    return fl_results

def getFlOnPoint(PyFLT flt_obj, float R, float Z, float Theta,
                 bool with_flt, int direction):
    """Obtain FL points using an origin point R, Z, Theta. Tracing in both
    directions.
    """

    cdef :
        vector[double] fl_buffer

    fl_buffer.clear()
    # Prepare the thread containers
    flt_obj.c_prepareThreadContainers()

    # Set initial parameters
    flt_obj.c_setIV(R, Z, Theta)
    # Set upward direction.
    flt_obj.c_setDirection(direction)

    flt_obj.c_getFL(fl_buffer, with_flt)

    # Do not cast <object> on fl_buffer. Cython automatically casts and
    # transforms the object to a python object.
    return fl_buffer
