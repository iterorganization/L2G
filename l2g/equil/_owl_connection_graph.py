import numpy as np
from math import isclose

class Polygon(object):
    """Polygon class for simple 2D Ray-Casting usage.
    """

    def __init__(self, poly_x: list[float], poly_y: list[float]):

        POLYX, POLYY = poly_x, poly_y

        N_POLYX = len(POLYX)
        # Remove any duplicate points.
        to_remove: list[int] = []
        for i in range(N_POLYX):
            next_i = (i+1) % N_POLYX
            if isclose(POLYX[i], POLYX[next_i]) and isclose(POLYY[i], POLYY[next_i]):
                to_remove.append(i)

        # Remove elements in reverse
        for el in to_remove[::-1]:
            POLYX.pop(el)
            POLYY.pop(el)

        self.poly_x: list[float] = POLYX
        self.poly_y: list[float] = POLYY
        self.N = len(self.poly_x)

    def checkIfOnEdge(self, px1: float, py1: float, px2: float, py2: float,
                      tx: float, ty: float) -> bool:
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
        # if isclose(px1, px2) and isclose(py1, py2):
        #     if isclose(tx, px1) and isclose(ty, py1):
        #         return True
        #     return False

        dx1 = tx - px1
        dy1 = ty - py1

        dx2 = px2 - px1
        dy2 = py2 - py1

        if isclose(dx1 * dy2 - dx2 * dy1, 0):
            c = True

        return c

    def checkIfIn(self, tx: float, ty: float, check_if_on_edge: bool = True) -> bool:
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

        N = self.N
        c = False
        j = N - 1
        for i in range(N):
            px1 = self.poly_x[i]
            py1 = self.poly_y[i]
            px2 = self.poly_x[j]
            py2 = self.poly_y[j]
            # First check if point is on the edge of these two points.
            if check_if_on_edge and self.checkIfOnEdge(px1, py1, px2, py2, tx, ty):
                return False

            # Use the recipe for Ray-Casting method, taken from:
            # https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
            if (py1 > ty) != (py2 > ty) and \
               (tx < (px2 - px1) * (ty - py1) / (py2 - py1) + px1):
               c = not c
            j = i
        return c

    def horizontalRayIntersection(self, tx: float, ty: float) -> float:
        """Function used to find where a horizontal line starting at point
        (tx, ty) hits the polygon. The polygon should encompass the point, if
        it does not this function will not tell you that but it will return
        the X position of the point.


        Arguments:
            tx (float): X position of the point
            ty (float): Y position of the point

        Returns:
            out (float): X position where the polygon is intersected. out is
                always > tx if the polygon encompasses the point (tx, ty).
        """
        out = tx
        N = self.N
        j = N - 1

        for i in range(N):
            px1 = self.poly_x[i]
            py1 = self.poly_y[i]
            px2 = self.poly_x[j]
            py2 = self.poly_y[j]

            if py2 > py1:
                test_y1, test_y2 = py1, py2
            else:
                test_y1, test_y2 = py2, py1

            if px2 < px1:
                px1, px2, py1, py2 = px2, px1, py2, py1

            if not (test_y1 <= ty <= test_y2):
                j = i
                continue
            out = px1 + (ty - py1) * (px2 - px1) / (py2 - py1)
            if out < tx:
                j = i
                continue
            else:
                break
        return out

def getOwlConlensGraph(eq: 'EQA'):
    """Use this function that obtains the OWL connection length of the input

    Arguments:
        eq (l2g.externals.equilibrium_analysis.EQA): EQA object containing an
            object with analyzed equilibrium data.

    Returns:
        drsep (np.ndarray): 1D array of distances from the boundary at the
            outer midplane
        conlen_down (np.ndarray): 1D array of connection length connecting
            outer midplane and outer divertor target
        conlen_up (np.ndarray): 1D array of connection length connecting
            outer midplane and inner divertor target
    """
    from l2g.external.bicubic import PyBicubic
    from l2g.external.rkf45 import PyRKF45FLT
    from math import sqrt, cos

    equilibrium = eq.equilibrium

    # Evaluate if still not evaluated
    eq.evaluate()

    Rb, Zc, Btotal, Bpm = eq.getMidplaneInfo(which="owl")

    # Get the R points on the midplane. Up to the wall.
    polygon = Polygon(equilibrium.wall_contour_r, equilibrium.wall_contour_z)
    r_end = polygon.horizontalRayIntersection(Rb, Zc)
    r_points = np.linspace(Rb, polygon.horizontalRayIntersection(Rb, Zc), 200)

    bicubic = PyBicubic(equilibrium.grid_r, equilibrium.grid_z,
                        equilibrium.psi)
    rkf45_obj = PyRKF45FLT()
    rkf45_obj.setInterpolator(bicubic)
    rkf45_obj.set_vacuum_fpol(equilibrium.fpol_vacuum)

    conlen_down = np.zeros(len(r_points))
    conlen_up = np.zeros(len(r_points))
    for i in range(len(r_points)):
        length = 0.0
        sr = r_points[i]
        sz = Zc
        sth = 0.0
        drsep = sr - Rb

        while length <= 500.0:
            nr, nz, nth = rkf45_obj.run_step(sr, sz, sth, 0.005)
            if not polygon.checkIfIn(nr, nz, False):
                break

            curr_l = sqrt(nr*nr + sr*sr - 2*nr*sr*cos(nth-sth) + (nz - sz)*(nz - sz))
            length += curr_l
            sr = nr
            sz = nz
            sth = nth
        conlen_down[i] = length


        # Now do it also in the opposite direction
        sr = r_points[i]
        sz = Zc
        sth = 0.0
        length = 0.0
        while length <= 500.0:
            nr, nz, nth = rkf45_obj.run_step(sr, sz, sth, -0.005)
            if not polygon.checkIfIn(nr, nz, False):
                break

            curr_l = sqrt(nr*nr + sr*sr - 2*nr*sr*cos(nth-sth) + (nz - sz)*(nz - sz))
            length += curr_l
            sr = nr
            sz = nz
            sth = nth
        conlen_up[i] = length
    return r_points - Rb, conlen_down, conlen_up

if __name__ == "__main__":
    import l2g
    import l2g.equil
    from l2g.external.equilibrium_analysis import EQA
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description="Python script for " +
        "obtaining the outer wall connection length by using the EQDSK-G " +
        "equilibrium wall contour for shadowing. It utilizes the external " +
        "C++ library.")

    parser.add_argument("eqdsk_file", help="EQDSK-G file")

    args = parser.parse_args()
    file = args.eqdsk_file

    if not os.path.exists(file):
        print(f"ERROR: File {file} does not exist!")
        sys.exit(1)

    equilibrium = l2g.equil.getEquilibriumFromEQFile("g900003_00230_ITER_15MA_eqdsk16HR.txt")
    import l2g.external.bfgs_2d
    l2g.enableDebugging()
    l2g.addStreamHandler()
    import logging
    log = logging.getLogger(__name__)

    eqa = EQA(equilibrium)
    eqa.evaluate()
    # eqa.getType()
    drsep, conlen_down, conlen_up = getOwlConlensGraph(eqa)

    import matplotlib.pyplot as plt

    plt.plot(drsep, conlen_down, label="down")
    plt.plot(drsep, conlen_up, label="up")
    plt.ylabel("Connection length [m]")
    plt.xlabel(r"$\Delta_{sep}$ - radial distane along the midplane [m]")
    plt.plot(drsep, conlen_down + conlen_up)
    plt.grid()
    plt.legend()
    plt.show()