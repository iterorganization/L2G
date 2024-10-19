from math import isclose

class Polygon(object):
    """Polygon class for checking if points are inside the polygon, or on the
    edge.
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

        If check_if_on_edge is True, the function becomes more time consuming.

        Arguments:
            polyX (list): X values of polygon vertices
            polyY (list): Y values of polygon vertices
            tx (float): X position of test point
            ty (float): Y position of test point
            check_if_on_edge (bool): Switch to also check if a point lies on
                edge.

        Returns:
            c (bool): Boolean, True, if point lies inside or on the edges of
                polygon else False.
        """

        N = self.N
        c = False
        j = N - 1
        for i in range(N):
            # First check if point is on the edge of these two points.
            px1 = self.poly_x[i]
            py1 = self.poly_y[i]
            px2 = self.poly_x[j]
            py2 = self.poly_y[j]
            if check_if_on_edge and self.checkIfOnEdge(px1, py1, px2, py2, tx, ty):
                return False

            # Use the recipe for Ray-Casting method, taken from:
            # https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
            if (py1 > ty) != (py2 > ty) and \
               (tx < (px2 - px1) * (ty - py1) / (py2 - py1) + px1):
               c = not c
            j = i
        return c

    def area(self) -> float:
        """Calculates the area of the (simple) polygon using the shoelace
        formula.
        """

        area = 0.0

        N = len(self.poly_x)

        for i in range(N):
            area += self.poly_x[i] * (self.poly_y[(i+1)%N] -
                                      self.poly_y[(i-1)%N])
        area *= 0.5
        return area