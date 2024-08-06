# distutils: language = c++
# cython: language_level = 3

# External CPP class
from l2g.external.flt cimport FLT
from l2g.external.embree cimport PyEmbreeAccell

from libcpp cimport bool
from libcpp.vector cimport vector

import logging
log = logging.getLogger()

cdef class PyFLT:
    """Python wrapper around the external FLT class. Used only for testing as
    the main FieldLineTracer class utilizes the external FLT class the most.
    """

    def __cinit__(self):
        self.c_flt = new FLT()

    def __init__(self):
        pass

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

    def getVacuumFPOL(self):
        """Returns the FPol in Vacuum.
        Returns:
            vacuumFpol (double): FPol in Vacuum. In m T.
        """
        return self.c_flt.getVacuumFPOL()

    def setShift(self, double rmove, double zmove):
        """Sets the radial and vertical shift of the plasma.

        Arguments:
            rmove (double): Radial shift in meters.
            zmove (double): Vertical shift in meters.
        """
        self.c_flt.setShift(rmove, zmove)

    def setAbsError(self, double abserr):
        """Sets the absolute error. Less than 1. Default 1e-5

        Arguments:
            abserr (double): Absolute error used by RKF45.
        """
        self.c_flt.setAbsError(abserr)

    def setRelError(self, double relerr):
        """Sets the relative error. Less than 1. Default 1e-5

        Arguments:
            abserr (double): relative error used by RKF45.
        """
        self.c_flt.setRelError(relerr)

    def setDesiredStep(self, double step):
        """Sets the parametric time or toroidal angle end and the resolution
        at which we wish gather points or to check for intersections of a FL.

        Arguments:
            step (double): Toroidal angle step at which we wish to gather
                           points. Used both for obtaining FL points and for
                           FLT intersection.
        """
        self.c_flt.setDesiredStep(step)

    def setMaximumFieldlineLength(self,double max_conlen):
        """Sets the maximum connection length to which we follow a FL.

        Arguments:
            max_conlen (double): Maximum connection length in meters.
        """
        self.c_flt.setMaximumFieldlineLength(max_conlen)

    def setSelfIntersectionAvoidanceLength(self,double value):
        """Sets the initial length at which we avoid intersection tests. For
        example 0.005 (meters), this means that at the first 5 millimeters, do
        not check for intersection tests.

        Arguments:
            value (double): Initial length at which intersection tests are
                            ignored. In meters.
        """
        self.c_flt.setSelfIntersectionAvoidanceLength(value)

    def prepare(self):
        """Calls the prepareInterpolation function of the FLT class.
        """
        cdef bool f
        f = self.c_flt.prepareInterpolation()
        return f

    def runFLT(self):
        """Calls the runFLT function which starts the FLT.
        """
        self.c_flt.runFLT()

    def getBCart(self, double r, double z, double phi):
        """Gets the Magnetic field vector in Cartesian coordinate system. Used
        when processing data on a geometry (getting quantities).

        Result is written in the out variable. (Bx, By, Bz)
        """
        cdef:
            vector[double] out
        self.c_flt.getBCart(r, z, phi, out)
        return <object> out

    def getBCyln(self, double r, double z):
        """Gets the Magnetic field vector (poloidal and toroidal component)
        in Cylindrical coordinate system. Used when processing data on a
        geometry (getting quantities).

        Result is written in the out variable. (Bpol, Btor)

        Arguments:
            r (double): Radial position
            z (double): Vertical position
            out (vector, inout): Vector for writing values
        """
        cdef:
            vector[double] out
        self.c_flt.getBCyln(r, z, out)
        return <object> out

    def getPoloidalFlux(self, double r, double z):
        """Gets the polodidal magnetic flux value. Used when processing data on
        a geometry (getting quantities).

        Returns:
            flux (double): Poloidal magnetic flux. In Webb/rad.

        Arguments:
            r (double): Radial position
            z (double): Vertical position
        """
        return self.c_flt.getPoloidalFlux(r, z)

    def applyRT(self, PyEmbreeAccell obj):
        """This function is used to bind an EmbreeAccell C++ object to the
        FLT C++ object. Only pointer is used, so essentially we could switch
        different Embree objects.
        """
        self.c_flt.setEmbreeObj(obj.c_eacc)
