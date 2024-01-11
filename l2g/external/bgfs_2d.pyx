# distutils: language = c++
# cython: language_level = 3

from l2g.external.bgfs_2d cimport PyBgfs2d
from l2g.external.bicubic cimport PyBicubic, BICUBIC_INTERP

# TODO:
# Stop when new position is outside the defined domain.


cdef class PyBgfs2d:
    def __cinit__(self):
        # Set the default bounds to something extremely large
        self.bx1 = -1e9
        self.bx2 = 1e9
        self.by1 = -1e9
        self.by2 = 1e9

    def __init__(self):
        pass

    def setInterpolator(self, PyBicubic obj):
        """Sets the Bicubic interpolator.
        """
        self.c_bicubic = obj.c_bicubic

    cdef void c_set_interpolator(self, BICUBIC_INTERP *interp):
        """Cython function for setting the bicubic object. Used only in other
        cython cdef code.
        """
        self.c_bicubic = interp

    def setBounds(self, x1: float, x2: float, y1: float, y2: float):
        # Set the bounds for the search area.
        self.bx1 = x1
        self.bx2 = x2
        self.by1 = y1
        self.by2 = y2

    cdef double line_search(self, double x, double y, double p1, double p2):
        """Backtrack line search with Wolfe conditions. In other words when
        trying to solve the problem of finding the minimum of an objective
        function, the line_search calculates an acceptable step length of a
        search direction that reduces the objective function sufficiently.

        The function is modified, that in the presence of a saddle point, it
        also moves into that direction and gets trapped.

        The conditions to determine a suitable step length:
         1.) The Armijo rule
            f(x_k + \alpha-k p_k) <= f(x_k) + c_1 \qlpha_k p_k^T \nabla f(x_k)
         2.) The curvature condition
            -p_k^T \nabla f(x_k + \alpha_k p_k) <= -c_2 p_k^T \nabla f(x_k)

        with 0 < c1 < c2 < 1, p_k is the descent direction. c_1 is usually
        chosen to be quite small while c_2 is larger.


        Arguments:
            x (float): Starting position x.
            y (float): Starting position y.
            p1 (float): Search direction in x.
            p2 (float): Search direction in y,

        """
        cdef:
            double a, c1, c2
            double fx, fdx, fdy
            double fx_new, fdx_new, fdy_new
            double x_new, y_new

            double dot1, dot2, dot3

        a = 1
        c1 = 1e-3
        c2 = 0.9

        self.c_bicubic.getValues(x, y, fx, fdx, fdy)

        # dot1 used for Armijo rule
        dot1 = fx + c1 * a * (fdx * p1 + fdy * p2)
        # dot3 is the rhs of the curvature condition
        dot3 = c2 * (fdx * p1 + fdy * p2)

        x_new = x + a * p1
        y_new = y + a * p2

        self.c_bicubic.getValues(x_new, y_new, fx_new, fdx_new, fdy_new)
        # dot 2 is the lhs of the curvature condition
        dot2 = fdx_new * p1 + fdy_new * p2

        #      Armijo rule      Curvature condition
        while fx_new >= dot1 or dot2 <= dot3:
            # Advance
            a = a * 0.5

            if a < 1e-4:
                # Hack to also get caught into saddle points. Which is then
                # checked with the second derivative test.
                a = 1
                break

            x_new = x + a * p1
            y_new = y + a * p2

            self.c_bicubic.getValues(x_new, y_new, fx_new, fdx_new, fdy_new)
            dot1 = fx + c1 * a * (fdx * p1 + fdy * p2)
            dot2 = fdx_new * p1 + fdy_new * p2

        return a

    def findMinimum(self, double guess_x, double guess_y, int max_it):
        """This is a gradient-descent minimizer function that  tries to find a
        point in the function domain space where the gradient is close to
        zero (the function also gets trapped into saddle points and instead of
        looping infinitely it returns the point as a extrema point).

        The main use of the function is finding the minimum of the function,
        hence the name remained, although it can get trapped and returns the
        position of a extrema point (saddle point).
        """
        cdef:
            double x, y
            double x_new, y_new
            double f, fdx, fdy
            double f_new, fdx_new, fdy_new
            double gradnorm # Gradient norm
            # Search direction coefficients (Newton method)
            double p1, p2
            # Hessian coefficients
            double h11, h12, h21, h22

            double lhr11, lhr12, lhr21, lhr22
            # Line search coefficient
            double a

            # Delta gradient coefficients
            double y1, y2
            # Movement coefficients
            double s1, s2
            # Leftover coefficients
            double r,
            double li11, li12, li21, li22
            double ri11, ri12, ri21, ri22
            Py_ssize_t it

            # Bounds values
            double bx1, bx2, by1, by2
            boold out_of_bounds

        # Bounds values
        bx1 = self.bx1
        bx2 = self.bx2
        by1 = self.by1
        by2 = self.by2
        out_of_bounds = False

        self.c_bicubic.getValues(guess_x, guess_y, f, fdx, fdy)

        # Initial Hessian
        h11 = 1.0
        h12 = 0.0
        h21 = 0.0
        h22 = 1.0

        gradnorm = fdx*fdx + fdy*fdy
        it = 0

        x = guess_x
        y = guess_y
        while gradnorm > 1e-10:
            if it > max_it:
                break
            it += 1

            # Direction vector. Opposite in the direction of the gradient.
            p1 = -(h11 * fdx + h12 * fdy)
            p2 = -(h21 * fdx + h22 * fdy)

            # Step estimation
            a = self.line_search(x, y, p1, p2)

            # Movement
            s1 = a * p1
            s2 = a * p2
            x_new = x + a * p1
            y_new = y + a * p2

            # Check if out of bounds
            if x_new < bx1 or x_new > bx2:
                out_of_bounds = True
                break
            if y_new < by1 or y_new > by2:
                out_of_bounds = True
                break

            self.c_bicubic.getValues(x_new, y_new, f_new, fdx_new, fdy_new)

            # Difference in gradients
            y1 = fdx_new - fdx
            y2 = fdy_new - fdy

            # Start calculating the update to the Hessian matrix.
            r = 1.0 / (y1*s1 + y2*s2)

            li11 = 1.0 - r * s1 * y1
            li12 = - r * s1 * y2
            li21 = - r * s2 * y1
            li22 = 1.0 - r * s2 * y2

            ri11 = 1.0 - r * y1 * s1
            ri12 = - r * y1 * s2
            ri21 = - r * y2 * s1
            ri22 = 1.0 - r * y2 * s2

            lhr11 = (li11 * h11 + li12 * h21) * ri11  + (li11 * h12 + li12 * h22) * ri21
            lhr12 = (li11 * h11 + li12 * h21) * ri12  + (li11 * h12 + li12 * h22) * ri22
            lhr21 = (li21 * h11 + li22 * h21) * ri11  + (li21 * h12 + li22 * h22) * ri21
            lhr22 = (li21 * h11 + li22 * h21) * ri12  + (li21 * h12 + li22 * h22) * ri22

            # Update Hessian matrix
            h11 = lhr11 + r * s1 * s1
            h12 = lhr12 + r * s1 * s2
            h21 = lhr21 + r * s2 * s1
            h22 = lhr22 + r * s2 * s2

            # Copy the new values.
            x = x_new
            y = y_new
            f = f_new
            fdx = fdx_new
            fdy = fdy_new
            # Calculate the norm of gradient
            gradnorm = fdx*fdx + fdy*fdy

        return out_of_bounds, x, y

    def secondDerivativeTest(self, double x, double y) -> bool:
        """Return the value of the second derivative test D.

        D = fxx fyy - fxy**2

        If D > 0 and fxx > 0:
            Local minimum at x
        If D > 0 fxx < 0:
            Local maximum at y
        If D < 0:
            Saddle point
        If D == 0:
            Unknown
        """

        cdef:
            double val, valdxdx, valdydy, valdxdy, dummy, D

        self.c_bicubic.getAllValues(x, y, val, dummy, dummy, valdxdy)
        self.c_bicubic.getSecondDerivatives(x, y, valdxdx, valdydy)

        return valdxdx * valdydy - valdxdy * valdxdy
