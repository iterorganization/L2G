from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize, bisect
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
import numpy as np

import logging
log = logging.getLogger(__name__)
from typing import Tuple

from l2g.equil import Equilibrium

class SplineInterpolator(object):
    """Wrapper around the RectBivariate spline interpolate algorithm.

    Reasoning:
        * Possibility of switching to a different interpolate algorithm.
        * The RectBivariateSpline takes it's argument (and when interpolating
            data) in the opposite order:
                - Interp[x, y, data] -> RectBivariateInterp[y, x, data]
                - Interp.evaluate(x_i, y_i) -> RectBivariateInterp.evaluate(y_i, x_i)
        * When dealing with the displacement, have it here, to avoid code
          duplication everywhere with the displacement offsets.
    """
    def __init__(self):
        self.Interp = None
        self.r_displ = 0.0
        self.z_displ = 0.0

    def interpolate_data(self, R, Z, data):
        """Calculates the splines over the data.
        """
        self.Interp = RectBivariateSpline(Z, R, data, kx=5, ky=5, s=0)

    def __call__(self, r, z, dr=0, dz=0, grid=False):
        return self.Interp(z - self.z_displ, r - self.r_displ, dx=dz, dy=dr, grid=grid)

    def ev(self, r, z, dr=0, dz=0):
        """Get the value of the interpolation at point (r, z).

        Arguments:
            r (float): R point. In meters.
            z (float): Z point. In meters.
            dr (float): Derivative order for R
            dz (float): Derivative order for Z
        """
        return self.Interp.ev(z - self.z_displ, r - self.r_displ, dx=dz, dy=dr)


def checkIfOnEdge(px1: float, py1: float, px2: float, py2: float, tx: float,
                   ty: float) -> bool:
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
    c = False

    # It's easier to remove points that are insanely close to each other prior
    # to calling this function
    # In the case the p1 and p2 is the same, ignore it.
    # if np.allclose(px1, px2) and np.allclose(py1, py2):
    #     if np.allclose(tx, px1) and np.allclose(ty, py1):
    #         return True
    #     return False

    dx1 = tx - px1
    dy1 = ty - py1

    dx2 = px2 - px1
    dy2 = py2 - py1

    if np.allclose(dx1 * dy2 - dx2 * dy1, 0):
        c = True

    return c

def checkIfPointInPoly(polyX: list, polyY: list, tx: float, ty: float) -> bool:
    """Checks if point (tx, ty) lies within or on the edges of polygon,
    consisting of vertices.

    Arguments:
        polyX (list): X values of polygon vertices
        polyY (list): Y values of polygon vertices
        tx (float): X position of test point
        ty (float): Y position of test point

    Returns:
        c (bool): Boolean, True, if point lies inside or on the edges of
            polygon else False.
    """

    N = len(polyX)
    c = False
    j = N - 1
    for i in range(N):
        # First check if point is on the edge of these two points.
        px1 = polyX[i]
        py1 = polyY[i]
        px2 = polyX[j]
        py2 = polyY[j]
        if checkIfOnEdge(px1, py1, px2, py2, tx, ty):
            return True

        # Use the recipe for Ray-Casting method, taken from:
        # https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
        if (py1 > ty) != (py2 > ty) and \
           (tx < (px2 - px1) * (ty - py1) / (py2 - py1) + px1):
           c = not c
        j = i
    return c


class EQ:
    """Class for evaluating values from a given equilibrium.
    """
    def __init__(self, equilibrium=None, r_displ=0.0, z_displ=0.0):
        self.resetValues()
        if equilibrium is not None:
            self.setEquilibrium(equilibrium)
        self.setDisplacement(r_displ, z_displ)

    def resetValues(self) -> None:
        # Plasma center
        # In EQDSK usually (RMAXIS, ZMAXIS)
        # ZMAXIS alse serves as the height of the artificially envisioned
        # midplane.

        # Have to check if one needs to get the "apex" of the Separatrix
        # curvature in order to determine the actual starting point of the
        # midplane
        self.oPoint = None

        # Plasma type: "div" or "lim"
        self.type_ = None

        # Limiter Point
        self.limPoint = None

        # X points
        self.lowXPoint = None
        self.uppXPoint = None

        # LCFS flux values.
        # If diverted, we can have two LCFS-es depending on what plasma profile
        # we use
        self.psiLCFS = None
        self.psiLCFS2 = None # Reserved for second separatrix

        self.contactPoint = None

        # From the LCFS flux values we can determine the boundary values, for
        # IWL or OWL:
        # Rm - radial position of the midplane point
        # B_total - magnetic magnitude at that point.
        # Bpm - poloidal component of the magnetic field.
        self._eq: Equilibrium = None

        self.r_displ = 0.0
        self.z_displ = 0.0

        self._verbose = 0

        self._psi_spline = SplineInterpolator()
        self._psi_center_sign = 1

        self.evaluated = False

        # Points of the LCFS contour, which are obtained when the type is
        # limiter.
        self._lcfs_points_1 = None
        self._lcfs_points_2 = None

    def setEquilibrium(self, obj: Equilibrium, load: bool = True,
                       reset: bool = True) -> None:
        """Set equilibrium data to the diagnostic class.

        Arguments:
            obj (Equilibirum): Equilibrium data (flux, fpol, ...)
            load (bool): Prepare interpolator and helper variables.
                Default True.
            reset (bool): Reset all the quantities used for analyzing the
                equilibrium. Useful when using single object for multiple
                equilibriums. Default True.

        """
        # First reset values. When a new equilibrium is added.
        if reset:
            self.resetValues()

        self._eq = obj
        self.evaluated = False
        if load:
            self.loadEqdskValues()

    def loadEqdskValues(self) -> None:
        if self._eq is None:
            return

        # Initial guess for Opoint
        self._initialGuessOPoint = [self._eq.mag_axis_r, self._eq.mag_axis_z]
        # self._OPointsBounds = [[3.5, 8.5], [-4.5, 4.7]]

        # Initial guess for Xpoint (lower)
        # TODO, set the initial guess position relative to the bounding box
        # of the wall silhouette.

        minR = np.min(self._eq.grid_r)
        maxR = np.max(self._eq.grid_r)
        minZ = np.min(self._eq.grid_z)
        maxZ = np.max(self._eq.grid_z)

        halfR = (maxR + minR) * 0.5
        lowerZ = minZ + 0.15 * (maxZ - minZ)
        upperZ = minZ + 0.85 * (maxZ - minZ)
        # Lower initial guess.
        self._initialGuessXPoint = [halfR ,lowerZ]

        # Initial guess for Xpoint (upper)
        # TODO, set the initial guess position relative to the bounding box
        # of the wall silhouette.
        self._initialGuessXPointUpper = [halfR, upperZ]

        self.interpolate()

        # Psi center sign used for scaning extremes where we use the minimize
        # function. Undo the displacement.
        self._psi_center_sign = np.sign(
            self._psi_spline.ev(self._eq.mag_axis_r, self._eq.mag_axis_z))

    def setDisplacement(self, r_displ, z_displ):
        self.r_displ = r_displ
        self.z_displ = z_displ
        self._psi_spline.r_displ = r_displ
        self._psi_spline.z_displ = z_displ

    def interpolate(self):
        self._psi_spline.interpolate_data(self._eq.grid_r, self._eq.grid_z,
                                          self._eq.psi)

        # self._psi_spline = RectBivariateSpline(self._eq.grid_z, self._eq.grid_r,
        #                                  self._eq.psi, kx=5, ky=5, s=0)

    def printValues(self, r, z):
        """Print values at Point (r, z) such as flux value, sum of square first
        partial derivatives and the secondary partial derivative test for
        saddle/extrema points.
        """

        flux = self._psi_spline(r, z)
        drdz_sq = self._psi_spline(r, z, dz=1)**2 + \
                  self._psi_spline(r, z, dr=1)**2
        dr2 = self._psi_spline(r, z, dr=2)
        dz2 = self._psi_spline(r, z, dz=2)
        drz = self._psi_spline(r, z, dr=1, dz=1)
        D = dr2 * dz2 - drz * drz
        log.info(f'Flux value at ({r}, {z}): {flux}')
        log.info(f'Sum of square first partial derivatives: {drdz_sq}')
        log.info(f'Second partial derivative test: {D}.')
        if D > 0:
            log.info('Local maximum')
        else:
            log.info('Saddle point')

    def scanExtrema(self, initialGuess = None):
        if not initialGuess:
            initialGuess = self._initialGuessOPoint
        sign = -1 if self._psi_center_sign > 0 else 1
        def evalPsi(p, *args):
            psi_spl = args[0]
            sign = args[1]
            return sign * psi_spl.ev(p[0], p[1])

        oPoint = minimize(evalPsi, initialGuess, args=(self._psi_spline, sign))
        # log.info(oPoint.message)
        return oPoint

    def scanSaddle(self, initialGuess = None):
        """Try to find a saddle point by searching the minimum of the sum of
        squares of first derivatives in the (R, Z) cylindrical plane.
        """

        if not initialGuess:
            initialGuess = self._initialGuessXPoint
        def evalPsiDxDy(p, *args):
            psi_spl = args[0]

            psi_dr1 = psi_spl.ev(p[0], p[1], dr=1)
            psi_dz1 = psi_spl.ev(p[0], p[1], dz=1)

            drdz = psi_dr1 * psi_dr1 + psi_dz1 * psi_dz1
            return drdz

        xPoint = minimize(evalPsiDxDy, initialGuess, args=(self._psi_spline))
        # log.info(xPoint.message)
        return xPoint

    def secondDerivativeTest(self, p):
        """Second derivative test to see if the point is a saddle. First the
        point is found by using the minimize function than this second test
        checks whether it is a saddle point.
        """
        # drdz_sq = self._psi_spline.ev(p[0], p[1], dr=1) **2 + \
        #     self._psi_spline.ev(p[0], p[1], dz=1) ** 2
        dr2 = self._psi_spline.ev(p[0], p[1], dr=2)
        dz2 = self._psi_spline.ev(p[0], p[1], dz=2)
        drz = self._psi_spline.ev(p[0], p[1], dr=1, dz=1)
        D = dr2 * dz2 - drz * drz

        return D > 0

    def getContourPaths(self, point=None, flux=None):
        f = plt.figure()
        ax = f.add_subplot(111)

        if flux is None:
            cs = ax.contour(self._eq.grid_r + self.r_displ,
                            self._eq.grid_z + self.z_displ,
                            self._eq.psi,
                            levels=[self._psi_spline(point[0], point[1])])
        else:
            cs = ax.contour(self._eq.grid_r + self.r_displ,
                            self._eq.grid_z + self.z_displ,
                            self._eq.psi,
                            levels=[flux])
        plt.close(f)

        paths = cs.collections[0].get_paths()
        outX = []
        outY = []
        for i in range(len(paths)):
            p = paths[i]
            v = p.vertices
            outX.append(v[:,0])
            outY.append(v[:,1])
        # del cs
        return outX, outY

    def get_midplane_info(self, lcfs:float=None, which:str='owl'):
        """Returns the information about the midplane. The midplane is
        evaluated by finding the extrema of the boundary contour.

        .. note::

           This function evaluates the boundary values, by first finding the
           height of the midplane by searching the radius extrema of the
           boundary contour.

        Calculates the Rb, Z, Btotal and Bpm.

        If no flux value is provided (boundary value) then the return result is
        a tuple of four -1 values.

        Arguments:
            lcfs (float): Boundary contour poloidal magnetic flux value
                (Web/rad)
            which (str): Which side: iwl (inner)| owl (outer). Default owl.


        Returns:
            (tuple): tuple containing:

              - Rb (float): Starting position in the major radius.
              - Z (float): Height of the midplane.
              - Btotal (float): Total magnetic field magnitude at start of the
                midplane
              - Bpm (float): Poloidal component at start of the midplane

        """

        if lcfs is None:
            lcfs = self.psiLCFS
        if lcfs is None:
            log.info('No flux value provided to determine boundary ' +
                         'values')
            return -1, -1, -1, -1


        # Initial position we start from on the boundary separatrix
        Z = self.oPoint[1] # Height of magnetic axis

        # Set the bisection borders.
        if which == 'iwl':
            a = self._eq.grid_r[0]
            b = self.oPoint[0]
        elif which == 'owl':
            a = self.oPoint[0]
            b = self._eq.grid_r[-1]
        else:
            log.info(f'Wrong side specified: {which}')
            return -1, -1, -1, -1

        def fun_bis(r, *args):
            spl, z, value = args
            return spl.ev(r, z) - value
        # Now find the interval so that we have valid border conditions for the
        # bisection method.

        r_points = np.linspace(a, b, 100)
        vals = [fun_bis(_, self._psi_spline, Z, lcfs) for _ in r_points]
        where_diff_sign = np.where(np.diff(np.sign(vals)) != 0)[0]

        if len(where_diff_sign) == 0:
            log.error("Cannot find boundary on the midplane!")
            return -1, -1, -1, -1

        if which == 'iwl':
            # Get the last hit
            a = r_points[where_diff_sign[-1]]
            b = r_points[where_diff_sign[-1] + 1]
        else: #owl
            # Get the first hit
            a = r_points[where_diff_sign[0]]
            b = r_points[where_diff_sign[0] + 1]

        # Now accurately get the radial position!
        Rb = bisect(fun_bis, a, b, args=(self._psi_spline, Z, lcfs))

        Bvec = self.getBCart(Rb, Z, 0)

        Btotal = np.linalg.norm(Bvec)

        # Poloidal component of the magnetic field
        psidr = self._psi_spline.ev(Rb, Z, dr=1)
        psidz = self._psi_spline.ev(Rb, Z, dz=1)
        Bpm = np.sqrt(psidr * psidr + psidz * psidz) / Rb

        # Toroidal component of the magnetic field.
        # Use the poloidal current function = F=Bt * R, where F is constant
        # outside the plasma core.
        Btor = self._eq.fpol_vacuum / Rb

        currentPsi = self._psi_spline.ev(Rb, Z)
        if not np.allclose(currentPsi, lcfs, 1e-5):
            log.warning("Extremal position on boundary contour is not tight with psiLCFS!")
            log.warning(f"{currentPsi} != {lcfs} (LCFS)")

        return Rb, Z, Btotal, Bpm

    def getOWL_midplane(self, lcfs=None) -> Tuple[float, float, float, float]:
        """See :py:meth:`get_midplane_info`
        """
        return self.get_midplane_info(lcfs=lcfs, which='owl')

    def getIWL_midplane(self, lcfs=None) -> Tuple[float, float, float, float]:
        """See :py:meth:`get_midplane_info`
        """
        return self.get_midplane_info(lcfs=lcfs, which='iwl')

    def getBCart(self, r, z, phi):
        """Returns the magnetic vector in Cartesian coordinate system.
        """
        # Get derivates.
        fluxdR = self._psi_spline.ev(r, z, dr=1)
        fluxdZ = self._psi_spline.ev(r, z, dz=1)

        br = - fluxdZ / r
        bphi = self._eq.fpol_vacuum / r
        bz = fluxdR / r
        return np.array([br * np.cos(phi) - bphi * np.sin(phi),
                br * np.sin(phi) + bphi * np.cos(phi),
                bz], dtype=np.float)

    def alignLcfsToPoint(self, point: np.ndarray) -> np.ndarray:
        """Finds the shortest displacement required to put the limiter LCFS
        to point p.

        Arguments:
            p (np.ndarray): 1D array with 2 values.

        Returns:
            displ (np.ndarray): New evaluated displacement.
        """

        if not self.type_ == "lim":
            log.info("alignLcfsToPoint not used for Diverted state")
            return

        if self._lcfs_points_1 is None or self._lcfs_points_2 is None:
            log.info("No LCFS points to use for aligning.")
            return

        diff_1 = point - self._lcfs_points_1
        dist_1 = np.sum(diff_1**2, axis=1)
        min_ind_1 = np.argmin(dist_1)

        diff_2 = point - self._lcfs_points_2
        dist_2 = np.sum(diff_2**2, axis=1)
        min_ind_2 = np.argmin(dist_2)

        # The new plasma displacement will point from the closest point on the
        # LCFS to the closest point on the geometry.
        if dist_1[min_ind_1] < dist_2[min_ind_2]:
            displ = diff_1[min_ind_1]
        else:
            displ = diff_2[min_ind_2]

        return displ

    def checkIfLimiter(self):
        """Checks if there is a contact point on the wall that is closed inside
        the wall.
        """
        # Limiter test, so we must find the highest Flux value on the wall
        # silhouette. Afterwards we follow a few points on this magnetic
        # surface to see if it is inside the wall silhouette
        log.info("Running EQ.checkIfLimiter()")
        Resolution = 0.001 # 1 millimeter


        # See if the Flux map is showing in the other direction
        wallMaxFlux = -1e10
        wallMaxPoint = None
        wallMinFlux = 1e10
        wallMinPoint = None


        RLIM, ZLIM = self._eq.wall_contour_r, self._eq.wall_contour_z

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

        # Apply the displacements too! But in this case we have to be
        # careful. As negative shift basically means that we shift the
        # plasma to the left. As well as the wall silhouette. As this is
        # the distance between the geometry and wall silhouette. Therefore
        # positive sign!!
        # RLIM = [el + self.r_displ for el in RLIM]
        # ZLIM = [el + self.z_displ for el in ZLIM]
        # Update. Do not add the displacements, otherwise they get applied
        # twice!!!!
        N = len(RLIM)

        RMAXIS, ZMAXIS = self._eq.mag_axis_r, self._eq.mag_axis_z
        c = 0

        log.info(f"Creating points on the limiter with resolution: {Resolution} m")
        for i in range(N - 1):
            divisor = RLIM[i+1] - RLIM[i]

            verticalFlag = 0
            if np.allclose(divisor, 0):
                verticalFlag = 1
            else:
                k = (ZLIM[i+1] - ZLIM[i]) / (RLIM[i+1] - RLIM[i])


            # Get points to analize

            if verticalFlag:
                Npoints = (ZLIM[i+1] - ZLIM[i]) / Resolution
                dR = Resolution
            else:
                n = ZLIM[i] - k * RLIM[i]
                dR = Resolution / np.sqrt(1 + k * k)
                Npoints = (RLIM[i + 1] - RLIM[i]) / dR

            if Npoints < 0:
                dR *= -1
            Npoints = int(abs(Npoints))

            for j in range(Npoints):
                c += 1
                if verticalFlag:
                    point = RLIM[i], ZLIM[i] + j * dR
                else:
                    point = RLIM[i] + j * dR, (RLIM[i] + j * dR) * k + n

                # if point[1] < -2.55:
                #     # Limit to avoid going into divertor area.
                #     continue
                flux = self._psi_spline.ev(point[0], point[1])

                # Check for monotone connection to Plasma center.
                # Tricky. If we move the plasma, of course the magnetic axis
                # moves with it! Therefore we need to move the RMAXIS and
                # ZMAXIS in the direction of the displacement.
                _arr = np.linspace((point[0], point[1]),
                                   (RMAXIS + self.r_displ, ZMAXIS + self.z_displ),
                                    20)
                _fluxArr = self._psi_spline(_arr[:, 0], _arr[:, 1])
                _fluxDiff = np.diff(_fluxArr)
                if np.all(_fluxDiff > 0.0) or np.all(_fluxDiff < 0.0):
                    # We have monotonic flux values to the center
                    if flux > wallMaxFlux:
                        wallMaxFlux = flux
                        wallMaxPoint = point
                    if flux < wallMinFlux:
                        wallMinFlux = flux
                        wallMinPoint = point
        if self._eq.psi_sign > 0:
            # The written boundary flux value is higher than in the center,
            # therefore we need to find the minimum flux value on the wall.
            contactFlux = wallMinFlux
            contactPoint = wallMinPoint
        else:
            # Otherwise find the max value
            contactFlux = wallMaxFlux
            contactPoint = wallMaxPoint

        log.info(f"Proposed contact value Psi [Wb] = {contactFlux}")
        log.info(f"Proposed contact position [m,m] = {contactPoint}")

        # print(wallMinFlux, wallMaxFlux)
        # print(contactFlux, contactPoint)
        # Now follow the curve somewhat and check if it is inside the chamber

        log.info("Tracking magnetic surface from contact point to see if it" +
                 " remains inside the wall.")
        def fun(t, p, *args):
            R, Z = p
            I_psi, sign = args
            valdr = -sign * I_psi.ev(R, Z, dz=1)
            valdz =  sign * I_psi.ev(R, Z, dr=1)
            norm = np.sqrt(valdr * valdr + valdz * valdz)

            return valdr / norm, valdz / norm
        sign = 1 # Follow in one of the direction
        solver = solve_ivp(fun=fun, t_span=(0, 2*np.pi), y0=contactPoint,
                           method='DOP853', t_eval=np.linspace(0,2*np.pi,100),
                           atol=1e-5, rtol=1e-5,
                           args=[self._psi_spline, sign])


        # import matplotlib.pyplot as plt
        # plt.plot(self._eq.wall_contour_r, self._eq.wall_contour_z)
        # plt.plot(solver.y[0], solver.y[1])
        # plt.axis('equal')
        # plt.show()
        self._lcfs_points_1 = solver.y.T

        # Experimental: find if the calculated trajectory make loops
        def get_loop_index(points: np.ndarray, atol: float=1e-3) -> int:
            index = -1
            diff_x = points[:, 0] - points[0, 0]
            diff_y = points[:, 1] - points[0, 1]
            diff = diff_x**2 + diff_y**2

            for i in range(1, len(diff)):
                if np.allclose(diff[i], 0.0, atol=atol):
                    index = i
                    break

            return index

        loop_index = get_loop_index(solver.y.T)
        if loop_index == -1:
            loop_index = len(solver.y[0])

        for i in range(loop_index):
            # print(solver.y[0, i], end='asd')
            t = checkIfPointInPoly(RLIM, ZLIM, solver.y[0, i], solver.y[1, i])
            if not t:
                log.info("Surface does not lie inside tokamak. Possibly " +
                         "diverted equilibrium.")
                return False, None, None

        # So for so good, let's save the lcfs points

        # Now the other direction. Sign changed at valx and valy
        sign = -1 # Now follow in the other direction
        solver = solve_ivp(fun=fun, t_span=(0, 2*np.pi), y0=contactPoint,
                           method='DOP853', t_eval=np.linspace(0,2*np.pi,100),
                           atol=1e-5, rtol=1e-5,
                           args=[self._psi_spline, sign])

        # import matplotlib.pyplot as plt
        # plt.plot(self._eq.wall_contour_r, self._eq.wall_contour_z)
        # plt.plot(solver.y[0], solver.y[1])
        # plt.axis('equal')
        # plt.show()
        log.info("Surface tracked, now to see if the points of the surface " +
                 "lies inside the tokamak.")
        # Ok now save the other points of the LCFS contour
        self._lcfs_points_2 = solver.y.T

        loop_index = get_loop_index(solver.y.T)
        if loop_index == -1:
            loop_index = len(solver.y[0])
        for i in range(loop_index):
            # print(solver.y[0, i], end='asd')
            t = checkIfPointInPoly(RLIM, ZLIM, solver.y[0, i], solver.y[1, i])
            if not t:
                log.info("Surface does not lie inside tokamak. Possibly " +
                         "diverted equilibrium.")
                return False, None, None


        log.info("Surface lies inside tokamak. Limiter configuration.")

        # If it is a limiter configuration, then the produced contour will make
        # two loops. Otherwise it is diverted!

        # Check bouncing in X axis.
        # xDiff = np.sign(np.diff(solver.y[0]))
        # print(xDiff)
        # c = 1
        # currentSign = xDiff[0]
        # for i in range(1, len(xDiff)):
        #     if xDiff[i] != currentSign:
        #         c += 1
        #         currentSign = xDiff[i]

        # if xDiff[0] == xDiff[-1]:
        #     c -= 1

        # print(f"Number of monotonic shifts: {c}")

        # if it is a limiter, the last point should be close to the starting
        # point
        # r0 = solver.y[0, 0] ** 2 + solver.y[1, 0] ** 2
        # r1 = solver.y[0, -1] ** 2 + solver.y[1, -1] ** 2
        # import matplotlib.pyplot as plt
        # plt.plot(self._eq.wall_contour_r, self._eq.wall_contour_z)
        # plt.plot(solver.y[0], solver.y[1])
        # plt.axis('equal')
        # plt.show()

        # if not np.allclose(r0, r1):
        #     return False, None, None

        # Now that we have a LCFS contour, let's save it for limiter cases, in
        # which we wish to position the LCFS directly on the target geometry,
        # by measuring the closest distance (x, y) from the closest point on
        # the target geometry to the LCFS.


        return True, contactFlux, contactPoint

    def calculate_drsep(self, sep1: float, sep2: float) -> float:
        """Calculates the drsep on the outer midplane, where the height of the
        midplane is the same as plasma center height.

        Only applicable for divertor type equilibriums.

        Returns:
            drsep (float): Distance between separatrixes, in mm
        """
        drsep = 0.0
        Z = self.oPoint[1] # Hight of magnetic axis

        # Use bisection to get radial positions of the separatrixes on the
        # midplane.


        def fun(r, *args):
            spl, z, value = args
            return spl.ev(r, z) - value

        # First manually check if drsep can be calculated.
        f1 = fun(self.oPoint[0], self._psi_spline, Z, sep1)
        f2 = fun(self._eq.grid_r[-1], self._psi_spline, Z, sep1)

        if f1 * f2 > 0.0:
            return 0.0

        f1 = fun(self.oPoint[0], self._psi_spline, Z, sep2)
        f2 = fun(self._eq.grid_r[-1], self._psi_spline, Z, sep2)

        if f1 * f2 > 0.0:
            return 0.0

        Rsep1 = bisect(fun, self.oPoint[0], self._eq.grid_r[-1],
                       args=(self._psi_spline, Z, sep1))
        Rsep2 = bisect(fun, self.oPoint[0], self._eq.grid_r[-1],
                       args=(self._psi_spline, Z, sep2))

        drsep = abs(Rsep2 - Rsep1)
        return 1e3 * drsep

    def evaluate(self):
        """Determine what type of equilibrium, get Psi fluxes.
        """
        # Check if it is limiter type or diverted
        if self.evaluated:
            log.info("Already evaluated.")
            return

        self.psiLCFS = None
        self.psiLCFS2 = None
        self.drsep = None

        log.info("Running EQ.evaluate()")
        log.info("First checking if equilibrium is limiter configuration.")
        isLim, contactFlux, contactPoint = self.checkIfLimiter()
        if isLim:
            self.type_ = "lim"
        else:
            self.type_ = "div"
        log.info(f"EQ.type_ = {self.type_}")

        # Search for O point.
        log.info("Trying to find O point.")
        oPoint = self.scanExtrema()

        # Check if derivatives are nearly zero

        dx1 = self._psi_spline.ev(oPoint.x[0], oPoint.x[1], dz=1)
        dy1 = self._psi_spline.ev(oPoint.x[0], oPoint.x[1], dr=1)
        dx2 = self._psi_spline.ev(oPoint.x[0], oPoint.x[1], dz=0)
        dy2 = self._psi_spline.ev(oPoint.x[0], oPoint.x[1], dr=2)
        v1 = dx1 * dx1 + dx2 * dx2
        v2 = dx2 * dx2 + dy2 * dy2
        if v1 < 1e-2 and v2 < 1e-2:
            pass
        else:
            log.debug('O point detected has too large derivatives:')
            log.debug(f'1st deriv sum square: {v1}')
            log.debug(f'2nd deriv sum square: {v2}')
        log.info(f'O point detected: {oPoint.x}')

        self.oPoint = oPoint.x

        if self.type_ == 'lim':
            self.psiLCFS = contactFlux
            self.contactPoint = contactPoint
            log.info('Limiter type detected:')
            log.info(f'Psi LCFS (Webb):{self.psiLCFS}')
            log.info(f'Contact point (R, Z) in meters: {self.contactPoint}')
            self.evaluated = True
            return

        log.info("Detecting X points")
        # Try to find Lower X-point, at the same time, this will be a
        # equilibrium type test
        lowX = self.scanSaddle()
        success = True
        if lowX.status:
            log.info("X search failed with non-zero status!")
            success = False
        else:
            log.info("X search succeeded, now performing 2nd derivative test")
            if self.secondDerivativeTest(lowX.x):
                log.info("2nd derivative test failed.")
                log.info(f"Detected X point is not a X point: {lowX.x}")
                log.info("No lower X point detected")
                success = False
            else:
                log.debug("2nd derivative test succeeded.")
                log.info("Obtained first separatrix flux value.")
                self.psiLCFS = self._psi_spline.ev(lowX.x[0], lowX.x[1])
                self.lowXPoint = lowX.x
                log.info(f"Psi 1 LCFS (Webb): {self.psiLCFS}")
                log.info(f"X point position (R, Z) in meters: {self.lowXPoint}")

        if not success:
            log.info("Failed to find X point in a diverted equilibrium!")
            log.info("Manually set the conditions or check if the input file is ok")

        # X point found, now let's use the 2nd derivative test. For
        # instance, it might go to the O point location

        log.info("Searching for secondary X point.")

        uppX = self.scanSaddle(self._initialGuessXPointUpper)

        if uppX.status:
            log.info("Upper X search failed with non-zero status")
        else:
            # X point found, now let's use the 2nd derivative test. For
            # instance, it might go to the O point location
            if self.secondDerivativeTest(uppX.x):
                # Local maximum
                log.debug("2nd derivative test failed.")
                log.info("Detected upper X point is not a saddle " +
                             f"point: {uppX.x}")
                self.drsep = None
            else:
                self.uppXPoint = uppX.x
                self.psiLCFS2 = self._psi_spline.ev(uppX.x[0], uppX.x[1])
                log.info(f"Psi 2 LCFS (Webb): {self.psiLCFS2}")
                log.info("X point position (R, Z) in meters: " +
                             f"{self.uppXPoint}")
                log.info("Both X points detected!")

        if self.psiLCFS and self.psiLCFS2:
            log.info("Calculating distance between the separatrixes.")
            self.drsep = self.calculate_drsep(self.psiLCFS, self.psiLCFS2)
            log.info(f"Drsep={self.drsep}")
        else:
            self.drsep = None

        # Now here is the tricky part. In some cases, we could have only one
        # upper X-point, like West...
        if self.psiLCFS is None and self.psiLCFS2:
            self.psiLCFS = self.psiLCFS2
            self.psiLCFS2 = None
            log.info("Upper X point detected, but lower not.")
            log.info("Assuming upper x point configuration.")


        self.evaluated = True
