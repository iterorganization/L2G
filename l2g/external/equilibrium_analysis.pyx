# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np
from libcpp cimport bool
from libc.stdio cimport printf
from libcpp.vector cimport vector
from libc.math cimport sqrt

import cython

from l2g.equil import Equilibrium

from l2g.external.bicubic cimport BICUBIC_INTERP
from l2g.external.rkf45 cimport RKF45
from l2g.external.bgfs_2d cimport PyBgfs2d

import logging
log = logging.getLogger(__name__)


cdef bool check_if_on_edge(const float px1, const float py1, const float px2,
                        const float py2, const float tx, const float ty):
    """Checks if point (tx, ty) lies on the edge. Using cross product instead
    of square root methods, which is much slower.

    Arguments:
        px1 (float) : X values of an edge vertice
        py1 (float) : Y values of an edge vertice
        px1 (float) : X values of an edge vertice
        py1 (float) : Y values of an edge vertice
        tx (float): X position of test point
        ty (float): Y position of test point

    Returns:
        c (bool): True if it lies on the edge, else False
    """
    cdef:
        bool c=False
        float dx1, dy1, dx2, dy2

    dx1 = tx - px1
    dy1 = ty - py1

    dx2 = px2 - px1
    dy2 = py2 - py1

    return abs(dx1 * dy2 - dx2 * dy1) < 1e-4

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef bool check_if_point_in_poly(const double[:] polyX, const double[:] polyY,
                                 const int polyN, const double tx,
                                 const double ty):

    cdef:
        bool c
        Py_ssize_t i, j
        double px1, px2, py1, py2
    c = False
    j = polyN - 1

    for i in range(polyN):
        px1 = polyX[i]
        py1 = polyY[i]

        px2 = polyX[j]
        py2 = polyY[j]

        if check_if_on_edge(px1, py1, px2, py2, tx, ty):
            return True
        # Use the recipe for Ray-Casting method, taken from:
        # https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
        if (py1 > ty) != (py2 > ty) and (tx < (px2 - px1) * (ty - py1) / (py2 - py1) + px1):
           c = not c
        j = i
    return c


cdef class EQA:
    cdef:
        BICUBIC_INTERP *c_bicubic
        RKF45 *c_rkf45
        PyBgfs2d c_bgfs
        bool evaluated
        object equilibrium
        str plasma_type
        tuple o_point, low_x_point, upp_x_point, contact_point
        double[:] wall_contour_r, wall_contour_z
        float psi_x, psi_upp_x, psi_lcfs, mag_axis_r, mag_axis_z, r_displ
        float min_z, max_z, min_r, max_r
        float z_displ, psi_grad_sign


    def __cinit__(self):
        self.c_bicubic = new BICUBIC_INTERP()
        self.c_bicubic.prepareContainers()

        self.c_rkf45 = new RKF45()
        self.c_rkf45.set_omp_thread(0)

    def __init__(self, equilibrium: Equilibrium):
        self.c_bgfs = PyBgfs2d()
        self.resetValues()
        self.setEquilibrium(equilibrium)

    def resetValues(self):
        self.evaluated = False
        self.equilibrium = None

        self.plasma_type = ""

        # Points
        self.o_point = None
        self.low_x_point = None
        self.upp_x_point = None
        self.contact_point = None

        # Flux characteristic value
        self.psi_x = 0.0
        self.psi_upp_x = 0.0
        self.psi_lcfs = 0.0
        self.wall_contour_r = None
        self.wall_contour_z = None
        self.mag_axis_r = 0.0
        self.mag_axis_z = 0.0

        # Displacement value
        self.r_displ = 0.0
        self.z_displ = 0.0

        # Direction of the gradient of psi
        self.psi_grad_sign = 0.0

    def setEquilibrium(self, obj: Equilibrium):
        """Sets the equilibrium objects and propagates the
        flux data to the interpolation and solver objects.
        """
        self.resetValues()

        self.equilibrium = obj
        self.evaluated = False

        # Set the arrays
        self.c_bicubic.setArrays(obj.grid_r,
                                 obj.grid_z, obj.psi)
        self.c_rkf45.set_interpolator(self.c_bicubic)
        self.c_rkf45.set_vacuum_fpol(self.equilibrium.fpol_vacuum)
        self.c_bgfs.c_set_interpolator(self.c_bicubic)

        # Clean the Wall contour
        RLIM, ZLIM = self.equilibrium.wall_contour_r, self.equilibrium.wall_contour_z

        N_RLIM = len(RLIM)
        # Remove any duplicate points.
        to_remove = []
        for i in range(N_RLIM):
            next_i = (i+1) % N_RLIM
            if np.allclose(RLIM[i], RLIM[next_i]) and np.allclose(ZLIM[i], ZLIM[next_i]):
                to_remove.append(i)

        # Remove elements in reverse
        for el in to_remove[::-1]:
            RLIM.pop(el)
            ZLIM.pop(el)

        self.wall_contour_r = np.array(RLIM)
        self.wall_contour_z = np.array(ZLIM)
        self.mag_axis_r = self.equilibrium.mag_axis_r
        self.mag_axis_z = self.equilibrium.mag_axis_z
        self.psi_grad_sign = self.equilibrium.psi_sign

    def setDisplacement(self, r_displ: double, z_displ: double):
        self.r_displ = r_displ
        self.z_displ = z_displ

    def getBoundaryFluxValue(self):
        if not self.plasma_type:
            log.error(f"No plasma type: {self.plasma_type}")
            return None

        if self.plasma_type == "div":
            return self.psi_x
        else:
            return self.psi_lcfs

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void get_contact_point(self, double *contact_flux, double *contact_r,
                                double *contact_z):
        """This function creates point on the wall silhouette and tries to
        find a contact point for a possible LCFS. The point has two conditions
        two fulfill, first it has to be the closest values (depending on the
        direction of the poloidal magnetic flux) and the values from the point
        to the magnetic axis have to be monotone (either strictly rising or
        falling).
        """
        cdef:
            double resolution = 0.001 # 1 millimeters
            double wall_max_flux = -1e10
            double wall_min_flux = 1e10

            double max_point_r, max_point_z
            double min_point_r, min_point_z

            double [:] RLIM, ZLIM
            double mag_axis_r, mag_axis_z

            double r_displ, z_displ

            # Various variables
            int i, j, k
            int N, n_points
            # Flag if a line is practically vertical.
            bool vertical_flag
            # Variables for creating points
            double divisor, slope, n
            double px, py

            double flux, dummy, prev_flux, next_flux

            # Monotone function check variables
            double mt_slope, mt_dx, mt_diff1, mt_diff2
            bool mt_monotone

        # Assign values to flux, dummy to squash warning
        flux = 0.0
        dummy = 0.0
        next_flux = 0.0

        r_displ = self.r_displ
        z_displ = self.z_displ
        # Read the magnetic axis position. Also take into account the
        # displacement values if any.

        RLIM = self.wall_contour_r
        ZLIM = self.wall_contour_z
        mag_axis_r = self.mag_axis_r + self.r_displ
        mag_axis_z = self.mag_axis_z + self.z_displ
        N = len(self.wall_contour_r)
        log.debug(f"Creating points on the limiter with resolution: {resolution} m")

        for i in range(N - 1):
            divisor = RLIM[i+1] - RLIM[i]
            vertical_flag = False
            if abs(divisor) < 1e-8:
                vertical_flag = True
            else:
                slope = (ZLIM[i+1] - ZLIM[i]) / divisor

            # Get points to analize
            if vertical_flag:
                n_points = int((ZLIM[i+1] - ZLIM[i]) / resolution)
                dr = resolution
                n = 1.0
            else:
                n = ZLIM[i] - slope * RLIM[i]
                # Same as cosine of the resolution length
                dr = resolution / sqrt(1 + slope * slope)
                n_points = int((RLIM[i + 1] - RLIM[i]) / dr)

            if n_points < 0:
                dr = -1.0 * dr

            n_points = int(abs(n_points))

            for j in range(n_points):
                if vertical_flag:
                    px = RLIM[i]
                    py = ZLIM[i] + j * dr
                else:
                    px = RLIM[i] + j * dr
                    py = (RLIM[i] + j * dr) * slope + n

                # Obtain flux values
                self.c_bicubic.getValues(px - r_displ, py - z_displ, flux,
                                         dummy, dummy)

                # Check for monotone connection to the plasma center.
                prev_flux = flux

                # Get the slope
                mt_slope = (mag_axis_z - py) / (mag_axis_r - px)
                mt_dx = (mag_axis_r - px) / 20.0
                self.c_bicubic.getValues(px + mt_dx - r_displ,
                                         py + mt_dx * mt_slope - z_displ,
                                         next_flux, dummy, dummy)
                mt_diff1 = next_flux - prev_flux
                prev_flux = next_flux
                mt_monotone = True
                for k in range(2, 21):
                    self.c_bicubic.getValues(px + k * mt_dx - r_displ,
                                             py + k * mt_dx * mt_slope - z_displ,
                                             next_flux, dummy, dummy)
                    mt_diff2 = next_flux - prev_flux

                    if mt_diff1 * mt_diff2 < 0.0:
                        mt_monotone = False
                        break

                if mt_monotone:
                    if flux > wall_max_flux:
                        wall_max_flux = flux
                        max_point_r = px
                        max_point_z = py
                    if flux < wall_min_flux:
                        wall_min_flux = flux
                        min_point_r = px
                        min_point_z = py

        # Essentially it would be good to make a check if there were any points
        # found
        if self.psi_grad_sign > 0:
            contact_flux[0] = wall_min_flux
            contact_r[0] = min_point_r
            contact_z[0] = min_point_z
        else:
            contact_flux[0] = wall_max_flux
            contact_r[0] = max_point_r
            contact_z[0] = max_point_z
        print("contact_flux %f\n", contact_flux[0])
        print("contact_r %f\n", contact_r[0])
        print("contact_z %f\n", contact_z[0])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void check_if_limiter(self, bool *is_closed, double *contact_r,
                               double *contact_z, double *contact_flux):
        log.debug("check_if_limiter()")

        cdef:
            # From equilibrium data memoeryview
            double [:] RLIM, ZLIM
            int N_RLIM

            # Rkf45 variables
            int flag
            double relerr, abserr, time, new_time, time_step
            double y[2]
            double yp[2]

        self.get_contact_point(contact_flux, contact_r, contact_z)
        print("contact_flux %f\n", contact_flux[0])
        print("contact_r %f\n", contact_r[0])
        print("contact_z %f\n", contact_z[0])
        if contact_r[0] == 0.0:
            log.debug("No contact point found")
            is_closed[0] = False
            return
        log.debug(f"Proposed contact value Psi [Wb/rad] = {contact_flux[0]}")
        log.debug(f"Proposed R,Z contact position [m,m] = {contact_r[0]} {contact_z[0]}")

        # Follow the magnetic surface to see if it is a limiter magnetic surface.
        relerr = 1e-4
        abserr = 1e-4

        # Do one loop around the tokamak.
        time_step = 2*np.pi/200

        y[0] = contact_r[0]
        y[1] = contact_z[0]
        time = 0.0
        flag = 1
        new_time = time_step

        RLIM = self.wall_contour_r
        ZLIM = self.wall_contour_z
        N_RLIM = len(RLIM)

        # Set displacements
        self.c_rkf45.set_r_move(self.r_displ)
        self.c_rkf45.set_z_move(self.z_displ)
        printf("I'm here: %f %f %f %f\n", y[0], y[1], time, new_time)
        for i in range(200):
            flag = self.c_rkf45.r8_rkf45(y, yp,  &time, new_time, &relerr, abserr, flag)
            printf("%f %f\n", y[0], y[1])
            new_time = time + time_step
            printf("asd")
            is_closed[0] = check_if_point_in_poly(RLIM, ZLIM, N_RLIM, y[0], y[1])
            printf("Seeing if check_point succeeded")
            if not is_closed[0]:
                break

    cpdef evaluate(self):
        """Evaluates the magnetic poloidal flux function.
        """
        cdef:
            bool is_closed=True
            double contact_r, contact_z, contact_flux

        # Use the information of the mag axis to find the O point.
        self.o_point = self.c_bgfs.findMinimum(self.equilibrium.mag_axis_r,
                                               self.equilibrium.mag_axis_z,
                                               100)

        contact_r = 0.0
        contact_z = 0.0
        contact_flux = 0.0


        if self.evaluated:
            log.debug("Already evaluated.")
            return
        if self.equilibrium is None:
            log.error("No equilibrium")

        self.psi_x = 0.0
        self.psi_upp_x = 0.0
        self.psi_lcfs = 0.0

        log.debug("Evaluating equilibrium")

        log.debug("Checking if it contains a closed flux surface")
        self.check_if_limiter(&is_closed, &contact_r, &contact_z, &contact_flux)
        if is_closed:
            log.debug("Closed flux surface found. Limiter configuration")
            self.plasma_type = "lim"
            self.contact_point = (contact_r, contact_z)
            log.debug("Finished evaluating")
            return

        # If we are here, it is not limiter!
        log.debug("No closed flux surface found. Diverted configuration")
        self.plasma_type = "div"

        # Now get the X-points!
        log.debug("Finding O point.")
        self.o_point = self.c_bgfs.findMinimum(self.equilibrium.mag_axis_r,
                                               self.equilibrium.mag_axis_z,
                                               100)
        log.debug(f"O point: {self.o_point}")
        minR = min(self.equilibrium.grid_r)
        maxR = max(self.equilibrium.grid_r)
        minZ = min(self.equilibrium.grid_z)
        maxZ = max(self.equilibrium.grid_z)
        half_r = (maxR + minR) * 0.5

        upper_z = minZ + 0.85 * (maxZ - minZ)
        self.upp_x_point = self.c_bgfs.findMinimum(half_r, upper_z, 100)
        log.debug(f"Upper X point: {self.upp_x_point}")

        # Setup the guesses for the X points
        lower_z = minZ + 0.15 * (maxZ - minZ)
        self.low_x_point = self.c_bgfs.findMinimum(half_r, lower_z, 100)
        log.debug(f"Lower X point: {self.low_x_point}")

        log.debug("Finished evaluating")

    def getType(self):
        return self.plasma_type

    def getContactPoint(self):
        return self.contact_point

    def getLowerXPoint(self):
        return self.low_x_point

    def getUpperXPoint(self):
        return self.upp_x_point
