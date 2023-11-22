# distutils: language = c++
# cython: language_level = 3

from l2g.external.bgfs_2d cimport PyBgfs2d
from l2g.external.bicubic cimport PyBicubic

cdef class PyBgfs2d:
    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def setInterpolator(self, PyBicubic obj):
        """Sets the Bicubic interpolator.
        """
        self.c_bicubic = obj.c_bicubic

    cdef double line_search(self, double x, double y, double p1, double p2):
        """Backtrack line search with Wolfe conditions. In other words

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
        dot1 = fx + c1 * a * (fdx * p1 + fdy * p2)
        dot3 = c2 * (fdx * p1 + fdy * p2)

        x_new = x + a * p1
        y_new = y + a * p2

        self.c_bicubic.getValues(x_new, y_new, fx_new, fdx_new, fdy_new)
        dot2 = fdx_new * p1 + fdy_new * p2
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
        """This function tries to find a point in the function domain space
        where the gradient is close to zero (the function also gets trapped
        into saddle points and instead of looping infinitely it returns the
        point as a extrema point).

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

            a = self.line_search(x, y, p1, p2)
            s1 = a * p1
            s2 = a * p2
            x_new = x + a * p1
            y_new = y + a * p2
            self.c_bicubic.getValues(x_new, y_new, f_new, fdx_new, fdy_new)

            y1 = fdx_new - fdx
            y2 = fdy_new - fdy

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

        return x, y

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
