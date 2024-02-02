# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.math cimport sqrt

import cython

from l2g.equil import Equilibrium

from l2g.external.bicubic cimport BICUBIC_INTERP
from l2g.external.rkf45 cimport RKF45
from l2g.external.bgfs_2d cimport PyBgfs2d

import logging
log = logging.getLogger(__name__)

@cython.wraparound(False)
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
        BICUBIC_INTERP *c_bicubic # For normal interpolation
        BICUBIC_INTERP *c_saddle_bicubic # For finding saddle points
        RKF45 *c_rkf45 # For tracing a magnetic surface
        PyBgfs2d c_bgfs
        bool evaluated
        object equilibrium
        str plasma_type
        list o_point, low_x_point, upp_x_point, contact_point, lcfs_points
        double[:] wall_contour_r, wall_contour_z
        float psi_x, psi_upp_x, psi_lcfs, mag_axis_r, mag_axis_z, r_displ
        float z_displ, psi_grad_sign
        object psi_dxdy


    def __cinit__(self):
        self.c_bicubic = new BICUBIC_INTERP()
        self.c_bicubic.prepareContainers()
        self.c_saddle_bicubic = new BICUBIC_INTERP()
        self.c_saddle_bicubic.prepareContainers()

        self.c_rkf45 = new RKF45()
        self.c_rkf45.set_omp_thread(0)

    def __init__(self, equilibrium: Equilibrium):
        self.c_bgfs = PyBgfs2d()
        self.resetValues()
        self.setEquilibrium(equilibrium)
        self.lcfs_points = []

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

        # Take the finite differences values array and use it for the saddle
        # points.
        # Remember: Z is the rows and R is the column
        self.psi_dxdy = np.zeros((obj.grid_z.size, obj.grid_r.size), dtype=float)
        for i in range(obj.grid_z.size):
            for j in range(obj.grid_r.size):
                self.psi_dxdy[i, j] = self.c_bicubic.m_fdx[i][j]*self.c_bicubic.m_fdx[i][j] + self.c_bicubic.m_fdy[i][j]*self.c_bicubic.m_fdy[i][j]
        self.c_saddle_bicubic.setArrays(obj.grid_r, obj.grid_z, self.psi_dxdy)
        self.c_rkf45.set_interpolator(self.c_bicubic)
        self.c_rkf45.set_vacuum_fpol(self.equilibrium.fpol_vacuum)
        self.c_bgfs.c_set_interpolator(self.c_bicubic)

        self.c_bgfs.setBounds(min(self.equilibrium.grid_r),
                              max(self.equilibrium.grid_r),
                              min(self.equilibrium.grid_z),
                              max(self.equilibrium.grid_z))

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
            int flag, i
            double relerr, abserr, time, new_time, time_step
            double y[2]
            double yp[2]
            double prev_distance

        self.get_contact_point(contact_flux, contact_r, contact_z)
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
        time_step = 0.01

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

        i = 0
        prev_distance = 0
        log.debug("Tracking mgnetic surface.")
        while i < 100000:
            flag = self.c_rkf45.r8_rkf45(y, yp,  &time, new_time, &relerr, abserr, flag)
            if flag != 2:
                # Redoing the step
                flag = 2
                continue
            new_time = time + time_step
            flag = 2
            i += 1
            # Store the points

            self.lcfs_points.append([y[0], y[1]])

            # Check if the point goes out of bounds.


            is_closed[0] = check_if_point_in_poly(RLIM, ZLIM, N_RLIM, y[0], y[1])
            if not is_closed[0]:
                log.debug("Magnetic surface not closed!")
                break

            # Get the distance from start
            if i > 100:
                # Also make sure that the distance is going down not above.
                d = (contact_r[0] - y[0])**2 + (contact_z[0] - y[1])**2

                if d < 0.1 and d < prev_distance: # less than 0.1 m^2, which corresponds to 0.01m
                    log.debug("Periodicity found after more than 100 steps! Stopped tracking magnetic surface")
                    break
                prev_distance = d
        log.debug(f"Magnetic surface followed. Is it closed? {is_closed[0]}")

    cpdef evaluate(self):
        """Evaluates the magnetic poloidal flux function.
        """
        cdef:
            bool is_closed=True
            double contact_r, contact_z, contact_flux

            double flux, dummy

        # Initializing flux and dummy to squash warnings
        flux = 0
        dummy = 0

        # Use the information of the mag axis to find the O point.
        out_of_bounds, *point = self.c_bgfs.findMinimum(self.equilibrium.mag_axis_r,
                                               self.equilibrium.mag_axis_z,
                                               100)
        if not out_of_bounds:
            self.o_point = point
            log.debug(f"O point: {self.o_point}")
        else:
            log.debug("Could not find O point!")

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
            self.psi_lcfs = contact_flux
            self.contact_point = [contact_r, contact_z]
            log.debug("Finished evaluating")
            self.evaluated = True
            return

        # If we are here, it is not limiter!
        log.debug("No closed flux surface found. Diverted configuration")
        self.plasma_type = "div"

        # Set the flux value interpolator
        self.c_bgfs.c_set_interpolator(self.c_bicubic)

        minR = min(self.equilibrium.grid_r)
        maxR = max(self.equilibrium.grid_r)
        minZ = min(self.equilibrium.grid_z)
        maxZ = max(self.equilibrium.grid_z)
        half_r = (maxR + minR) * 0.5

        # Set the derivative points interpoaltor
        self.c_bgfs.c_set_interpolator(self.c_saddle_bicubic)

        upper_z = minZ + 0.85 * (maxZ - minZ)
        log.debug("")
        log.debug("Upper X point search")
        out_of_bounds, *point = self.c_bgfs.findMinimum(half_r, upper_z, 1000)
        if not out_of_bounds:
            self.upp_x_point = point
            self.c_bicubic.getValues(self.upp_x_point[0] - self.r_displ,
                                     self.upp_x_point[1] - self.z_displ,
                                     flux,
                                     dummy, dummy)
            self.psi_upp_x = flux
            log.debug(f"Upper X point: {self.upp_x_point}")
        else:
            log.debug(f"Could not find upper X point!")

        # Setup the guesses for the X points
        lower_z = minZ + 0.15 * (maxZ - minZ)
        log.debug("")
        log.debug("Lower X point search")
        out_of_bounds, *point = self.c_bgfs.findMinimum(half_r, lower_z, 1000)

        if not out_of_bounds:
            self.low_x_point = point
            log.debug(f"Lower X point: {self.low_x_point}")
            self.c_bicubic.getValues(self.low_x_point[0] - self.r_displ,
                                     self.low_x_point[1] - self.z_displ,
                                     flux,
                                     dummy, dummy)
            self.psi_x = flux
        else:
            log.debug(f"Could not find lower X point!")
        # This is also the location of the separatrix
        log.debug("Finished evaluating")
        self.evaluated = True

    def getType(self):
        return self.plasma_type

    def getBoundaryFluxValue(self):
        if not self.plasma_type:
            # log.error(f"No plasma type: {self.plasma_type}")
            return None

        if self.plasma_type == "div":
            return self.psi_x
        else:
            return self.psi_lcfs

    def getSecondaryXFluxValue(self):
        if not self.plasma_type:
            # log.error(f"No plasma type: {self.plasma_type}")
            return None
        if not self.plasma_type == "div":
            log.error(f"Plasma is not diverted type or does not have upper X point")
        return self.psi_upp_x

    def getContactPoint(self):
        return self.contact_point

    def getLowerXPoint(self):
        return self.low_x_point

    def getUpperXPoint(self):
        return self.upp_x_point

    cdef double bisection(self, double flux, double Z, double a, double b):
        cdef:
            double c, tol
            double vala, valb, valc, dummy
            int n, max_iter

        log.debug(f"Borders a={a}, b={b}, flux={flux}, Z={Z}")
        tol = 1e-6

        n = 0
        max_iter = 100

        dummy = 0
        vala = 0
        valb = 0
        valc = 0

        self.c_bicubic.getValues(a - self.r_displ, Z - self.z_displ,
                                 vala, dummy, dummy)
        self.c_bicubic.getValues(b - self.r_displ, Z - self.z_displ,
                                 valb, dummy, dummy)

        log.debug("Starting bisection method")
        while n < max_iter:
            # Get the midpoint
            c = 0.5 * (a + b)

            self.c_bicubic.getValues(c - self.r_displ, Z - self.z_displ,
                                     valc, dummy, dummy)

            if abs(valc - flux) == 0 or 0.5 * (b - a) < tol:
                log.debug(f"Bisection finished after {n} steps, c={c}")
                return c

            n += 1
            if (valc - flux) * (vala - flux) > 0: # Same sign
                a = c
                vala = valc
            else:
                b = c
                valb = valc

        log.error("Bisection failed.")
        return 0

    def getMidplaneInfo(self, lcfs: float=None, which:str='owl'):
        cdef:
            # For the bisection
            double a, b # For the borders
            double st, da # For finding the borders
            # For the magnetic component
            double psidr, psidz
            double br, bphi, bz, Btotal, Bpm
            double Rb, Z, flux
            bool found_borders
            # For the interpolation
            double dummy, fl, pr_fl
            Py_ssize_t i

        # Initialize the values to suppress the warnings
        flux = 0
        dummy = 0
        fl = 0
        pr_fl = 0
        psidr = 0
        psidz = 0
        found_borders = False

        if lcfs is None:
            if self.plasma_type == "div":
                flux = self.psi_x
            else:
                flux = self.psi_lcfs

        Z = self.o_point[1] # Height of the magnetic axis

        # Find the interval so that we have a valid border conditions for the
        # bisection method.

        if which == "iwl":
            a = self.equilibrium.grid_r[0]
            b = self.o_point[0]
            da = -(b - a) / 100.0
            st = b
        elif which == "owl":
            a = self.o_point[0]
            b = self.equilibrium.grid_r[-1]
            da = (b - a) / 100.0
            st = a
        else:
            log.error(f"Wrong side specified: {which}")
            return -1, -1, -1, -1


        self.c_bicubic.getValues(st - self.r_displ, Z - self.z_displ,
                                 pr_fl, dummy, dummy)

        for i in range(1, 101):
            self.c_bicubic.getValues(st + i * da - self.r_displ,
                                     Z - self.z_displ, fl, dummy, dummy)

            if (pr_fl - flux) * (fl - flux) < 0.0:
                b = a + i * da
                a = b - da
                found_borders = True
                break

        if not found_borders:
            log.error("Cannot find the magnetic surface on the midplane")
            return -1, -1, -1, -1
        # Now get the major radius of the magnetic surface
        Rb = self.bisection(flux, Z, a, b)

        self.c_bicubic.getValues(Rb - self.r_displ, Z - self.z_displ,
                                 dummy, psidr, psidz)


        # Calculate the poloidal component of the magnetic field
        Bpm = sqrt(psidr * psidr + psidz * psidz) / Rb
        # Not really needed, Toroidal component of the magnetic field
        # Btor = self.equilibrium.fpol_vacuum / Rb

        # Calculate the total component of the magnetic field.
        br = psidz / Rb
        bphi = self.equilibrium.fpol_vacuum / Rb
        bz = psidr / Rb
        Btotal = sqrt(br * br + bphi * bphi + bz * bz)

        return Rb, Z, Btotal, Bpm

    def alignLcfsToPoint(self, point: list) -> list:
        """Finds the shortest displacement required to put the limiter LCFS
        to point p.

        Arguments:
            p (list): 1D array with 2 values.
        Returns:
            displ (np.ndarray): New evaluated displacement.
        """

        cdef:
            int i, mini
            double dist, min_dist
            list displacement

        if not self.plasma_type == "lim":
            log.error("alignLcfsToPoint only used for Limiter state.")
            return

        if not self.evaluated:
            log.error("Equilibrium not yet evaluated. Run .evaluate() first")
            return

        if not self.lcfs_points:
            log.error("No LCFS points stored in order to find the shortest distance.")
            return
        mini = -1
        min_dist = 1e6
        for i in range(len(self.lcfs_points)):
            dist = (self.lcfs_points[i][0] - point[0]) ** 2 + \
                   (self.lcfs_points[i][1] - point[1])

            if dist < min_dist:
                min_dist = dist
                mini = i

        return [point[0] - self.lcfs_points[mini][0], point[1] - self.lcfs_points[mini][1]]

    def distanceBetweenPsiOnMidplane(self, sep1: float, sep2: float) -> float:
        """Calculates the difference between the provided magnetic surface on
        the midplane.

        Returns:
            drsep (float): Distance between separatrixes, in mm
        """

        cdef:
            double drsep, Z
            # For the bisection
            double a, b # For the borders
            double st, da # For finding the borders
            bool found_borders
            # For the magnetic component
            double fl, pr_fl, dummy

            double sep1r, sep2r

        # Variable initialization
        drsep = 0
        fl = 0
        pr_fl = 0
        dummy = 0
        Z = self.o_point[1]

        # Find the borders for for sep1 and then sep2
        a = self.o_point[0]
        b = self.equilibrium.grid_r[-1]
        st = a
        da = (b - a) / 100.0

        # Get the first distance
        self.c_bicubic.getValues(a - self.r_displ, Z - self.z_displ, pr_fl,
                                 dummy, dummy)
        found_borders = False
        for i in range(1, 101):
            self.c_bicubic.getValues(st + i * da - self.r_displ,
                                     Z - self.z_displ, fl, dummy, dummy)

            if (pr_fl - sep1) * (fl - sep1) < 0.0:
                b = a + i * da
                a = b - da
                found_borders = True
                break
        if not found_borders:
            log.error("Cannot find the magnetic surface on the midplane")
            return 0

        sep1r = self.bisection(sep1, Z, a, b)

        # Find the borders for for sep1 and then sep2
        a = self.o_point[0]
        b = self.equilibrium.grid_r[-1]
        st = a
        da = (b - a) / 100.0

        # Get the first distance
        self.c_bicubic.getValues(a - self.r_displ, Z - self.z_displ, pr_fl,
                                 dummy, dummy)
        found_borders = False
        for i in range(1, 101):
            self.c_bicubic.getValues(st + i * da - self.r_displ,
                                     Z - self.z_displ, fl, dummy, dummy)

            if (pr_fl - sep2) * (fl - sep2) < 0.0:
                b = a + i * da
                a = b - da
                found_borders = True
                break
        if not found_borders:
            log.error("Cannot find the magnetic surface on the midplane")
            return 0

        sep2r = self.bisection(sep2, Z, a, b)
        return abs(sep1r - sep2r) * 1e3

    def getPsiDxDy(self):
        return self.psi_dxdy