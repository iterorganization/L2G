# distutils: language = c++
# cython: language_level = 3

from l2g.external.bfgs_2d cimport PyBfgs2d
from l2g.external.bicubic cimport PyBicubic, BICUBIC_INTERP

from libc.math cimport sqrt

from libcpp cimport bool

import logging
log = logging.getLogger(__name__)

cdef class PyBfgs2d:
    def __cinit__(self):
        # Set the default bounds to something extremely large
        self.bx1 = -1e9
        self.bx2 = 1e9
        self.by1 = -1e9
        self.by2 = 1e9

    def __init__(self):
        self._store = []

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

    cdef double quadmin(self, double a, double fa, double fpa, double b, double fb):
        """Finds the minimum solution for a quadratic polynomial.


        f(x) = A (x - a)**2 + B (x - a) + C


        To get the xmin where f'(xmin) = 0.0, do the following:

        1.) C = fa, through f(a) = fa
        2.) B = fpa, through f'(a) = fpa
        3.) A = (fb - fp(b-a) - fa) / (b-a)**2

        Get xmin through f'(xmin) = 2A(x-a) + B = 0
        """

        cdef:
            double A

        A = (fb - fpa*(b-a) - fa) / ((b-a)*(b-a))

        if A == 0.0:
            return 0.0
        return a - fpa / (2.0*A)


    cdef double cubicmin(self, double a, double fa, double fpa, double b, double fb, double c, double fc):
        """Finds the minimum of the cubic polynomial of

        f(x) = A (x - a)**3 + B (x - a)**2 + C (x-a) + D

        # In order to get the x where the function value is minimal:

        1.) Get D from f(a) = fa
        2.) Get C from f'(a) = fpa
        3.) Get A, B from solving the linear system of equations using points
            (b, fb) and (c, fc)

        Now finally solve f'(xmin) = 0, by solving the quadratic formula

        3A(x-a)**2 + 2B(x-a) + C = 0
        """

        cdef:
            double lb, lc # Lengths between b|c and a.
            double denom

            # Only variables A and B have to be computed
            double A, B

            # Variables used to solve the 2x2 linear system of equations
            double d11, d12, d21, d22
            double c1, c2
            double D # From solving the quadratic equation

        # Calculate the A and B
        lb = b - a
        lc = c - a

        # Inverse the 2x2 matrix.
        denom = lb * lc * lb * lc * (lb - lc)
        if denom == 0.0:
            return 0.0 # Bisection takes over

        d11 = lc*lc
        d12 = -lb*lb
        d21 = -lc * lc * lc
        d22 = lb * lb * lb
        c1 = fb - fa - fpa * lb
        c2 = fc - fa - fpa * lc

        A = d11 * c1 + d12 * c2

        if A == 0.0:
            return 0.0 # Bisection takes over
        B = d21 * c1 + d22 * c2

        A /= denom
        B /= denom

        D = B * B - 3.0 * A * fpa
        if D < 0.0:
            return 0.0 # Bisection takes over

        return a + (-B + sqrt(B * B - 3.0 * A * fpa)) / (3.0 * A)


    cdef double zoom(self, double x, double y, double a1, double a2,
                     double f0, double gk0, double p1, double p2, double c1,
                     double c2):
        """Zoom function of the approximate line-search satisfying strong
        Wolfe conditions.

        The process is as following:

            We have a Phi function Phi(a) = f(x + ap1, y + a p2), for which
            we want to find the step length a, that minimizes the function
            adequately along a direction vector vec_p=(p1, p2). The step length
            must be such that it satisfies the strong Wolfe conditions (Armijo
            and curvature conditions).

            Given an interval [a1, a2] where the starting a1 is very low,
            (i.e., 0) which always satisfies the line-search conditions and
            a2 is a value which does not satisfies the line-search
            conditions, find the step-length a* that adequately minimizes the
            function but not overshoots too much when going into that
            direction.

            1.) Use one of the methods, quadratic, cubic minimizers
                or simply bisection, that try to find a a* which minimizes the
                function the most and see if it satisfies the strong Wolfe
                condition.

            2.) If the calculated a* fails to satisfies the (normal) Wolfe
                conditions, we zoom to the [a1, a*] interval. We go back to
                step 1.)

            3.) (Normal) Wolfe conditions are satisfied.
                a.) We check for strong Wolfe conditions. If they are satisfied
                    we return a* as the proposed step length.
                b.) We zoom to a new interval. In this case we look at the dot
                    product between the gradient of Phi(a*) and the
                    direction vector to decide the new interval. Go to step
                    1.).

        """

        cdef:
            int maxiter, i
            # Cubic and quadratic interpolant check
            double delta1, delta2
            double quadchk, cubichk, gkinterp
            double da, aj
            double ahi, alo, ac
            double f1, fdx1, fdy1
            double f2, fdx2, fdy2
            double fj, fdxj, fdyj
            double fc, dummy
            double gkj

        # Initialize variables to squash referenced before assignment warning
        f1 =0.0
        fdx1 = 0.0
        fdy1 = 0.0
        f2 = 0.0
        fdx2 = 0.0
        fdy2 = 0.0
        fj = 0.0
        fdxj = 0.0
        fdyj = 0.0
        fc = 0.0
        dummy = 0.0

        delta1 = 0.2
        delta2 = 0.1
        maxiter = 10
        ac = 0.0
        i = 0
        # The dot product between the direction vector and derivative.

        ahi = a2
        alo = a1
        ac = 0.0

        while 1:
            # Get value at borders
            self.c_bicubic.getValues(x + alo * p1, y + alo * p2, f1, fdx1, fdy1)
            self.c_bicubic.getValues(x + ahi * p1, y + ahi * p2, f2, fdx2, fdy2)
            self.c_bicubic.getValues(x +  ac * p1, y +  ac * p2, fc, dummy, dummy)

            da = ahi - alo
            # a2 and a1
            if da < 0.0:
                a2 = alo
                a1 = ahi
            else:
                a1 = alo
                a2 = ahi


            # Start with quadratic minimizing
            if i == 0:
                # Use the quadratic interpolation minimizer method using the
                # alpha_hi, alpha_lo values.
                quadchk = delta2 * da
                gkinterp = fdx1 * p1 + fdy1 * p2
                aj = self.quadmin(alo, f1, gkinterp, ahi, f2)
                if aj > a2 - quadchk or aj < a1 + quadchk:
                    # Too close to the borders or not successful
                    aj = alo + 0.5 * da
            else:
                # Use the cubic interpolation minimizer method using the
                # alpha_hi, alppha_lo and one more value obtained when
                # switching the interval.
                cubichk = delta1 * da
                gkinterp = fdx1 * p1 + fdy1 * p2
                aj = self.cubicmin(alo, f1, gkinterp, ahi, f2, ac, fc)
                if aj > a2 - cubichk or aj < a1 + cubichk:
                    # Too close to the borders or not successful
                    aj = alo + 0.5 * da
            # Let's look at aj now
            self.c_bicubic.getValues(x + aj * p1, y + aj * p2, fj, fdxj, fdyj)

            # Zoom in to lower interval if we fail to check conditions
            if fj > f0 + c1 * aj * gk0 or fj >= f1:
                ac = a2
                ahi = aj
            else:
                # See if strong Wolfe condition is accepted
                gkj = fdxj*p1 + fdyj*p2
                if abs(gkj) <= -c2 * gk0:
                    # We found our winner
                    return aj

                # Move the interval ends. ahi or alo values do not necessarily
                # correspond to the lower, upper values of the interval
                # respectively as you would normally expect. This is remedied
                # with the sign of da (distance between ahi and alo).
                if gkj * da >= 0.0:
                    # Save the values for the cubic interpolation
                    ac = ahi
                    ahi = alo
                else:
                    ac = alo
                alo = aj

            i += 1
            if i > maxiter:
                log.debug("Failed to find a conforming step size.")
                return aj
        return 0.0

    cdef double line_search(self, double x, double y, double p1, double p2):
        """Backtrack line search algorithm with Wolfe conditions and zooming
        with cubic and quadratic minimizers.

        The post fixes 1 and 2 in the variable names signifies left and right
        end of the search interval respectfully.

        See Wright and Nocedal, 'Numerical Optimization', 1999
        """

        cdef:
            double a1, a2, aout
            double f0, fdx0, fdy0
            double f1, fdx1, fdy1
            double f2, fdx2, fdy2
            double gk0, gk2
            double dummy
            double c1, c2
            int maxiter, i

        # Values at the starting x, y
        f0 = 0.0
        fdx0 = 0.0
        fdy0 = 0.0

        f1 = 0.0
        fdx1 = 0.0
        fdy1 = 0.0
        f2 = 0.0
        fdx2 = 0.0
        fdy2 = 0.0

        a1 = 0
        a2 = 1.0
        aout = 1.0
        c1 = 1e-4
        c2 = 0.9
        maxiter = 10

        # Get the function value at the start and end
        self.c_bicubic.getValues(x, y, f0, fdx0, fdy0)
        gk0 = fdx0 * p1 + fdy0 * p2
        f1 = f0
        fdx1 = fdx0
        fdy1 = fdy0

        self.c_bicubic.getValues(x + a2 * p1,
                                 y + a2 * p2,
                                 f2, fdx2, fdy2)

        for i in range(maxiter):
            if a2 == 0.0:
                # Usually this should not happen.
                log.error("Line search failed with upper bound in line search falling to zero")
                return 0.0

            # Dot product between the gradient and direction vector.
            gk2 = fdx2 * p1 + fdy2 * p2

            # Invoke the zooming by checking if Wolfe's condition for function
            # decrease is not holding or if the function value at the end
            # of the length interval is higher than at the beginning.
            if f2 > f0 + c1 * a2 * gk0 or (f2 >= f1 and i > 0):
                aout = self.zoom(x, y, a1, a2, f0, gk0, p1, p2, c1, c2)
                break


            # Strong Wolfe condition.
            if abs(gk2) <= -c2 * gk0:
                return a2

            # Invoke zooming if the dot product of the direction and the
            # derivative values on the end of the a interval is > 0
            if gk2 >= 0:
                aout =  self.zoom(x, y, a1, a2, f0, gk0, p1, p2, c1, c2)
                break

            # Increase interval ends.

            a1 = a2
            a2 = 2 * a2
            # Update F values.
            self.c_bicubic.getValues(x + a1 * p1, y + a1 * p2, f1, fdx1, fdy1)
            self.c_bicubic.getValues(x + a2 * p1, y + a2 * p2, f2, fdx2, fdy2)

        else:
            # Maximum iteration reached
            log.debug("Maximum iteration reached. Stopping and setting a=1.0")
            return 1.0
        return aout

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
            bool out_of_bounds

        # Clear stored points
        self._store = [[guess_x, guess_y]]

        # Bounds values
        bx1 = self.bx1
        bx2 = self.bx2
        by1 = self.by1
        by2 = self.by2
        out_of_bounds = False

        f = 0.0
        fdx = 0.0
        fdy = 0.0
        f_new = 0.0
        fdx_new = 0.0
        fdy_new = 0.0

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

        log.debug(f"Starting the minimization method from point {x},{y}")
        while gradnorm > 1e-10:
            if it> max_it:
                log.debug("Exceeded maximum iterations")
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
                log.debug("Went outside of defined X bounds")
                break
            if y_new < by1 or y_new > by2:
                log.debug("Went outside of defined Y bounds")
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

            # lhr11 = (li11 * h11 + li12 * h21) * ri11  + (li11 * h12 + li12 * h22) * ri21
            # lhr12 = (li11 * h11 + li12 * h21) * ri12  + (li11 * h12 + li12 * h22) * ri22
            # lhr21 = (li21 * h11 + li22 * h21) * ri11  + (li21 * h12 + li22 * h22) * ri21
            # lhr22 = (li21 * h11 + li22 * h21) * ri12  + (li21 * h12 + li22 * h22) * ri22

            lhr11 = li11 * (h11 * ri11 + h12 * ri21) + li12 * (h21 * ri11 + h22 * ri21)
            lhr12 = li11 * (h11 * ri12 + h12 * ri22) + li12 * (h21 * ri12 + h22 * ri22)
            lhr21 = li21 * (h11 * ri11 + h12 * ri21) + li22 * (h21 * ri11 + h22 * ri21)
            lhr22 = li21 * (h11 * ri12 + h12 * ri22) + li22 * (h21 * ri12 + h22 * ri22)
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
            self._store.append([x, y])
            # Calculate the norm of gradient
            gradnorm = fdx*fdx + fdy*fdy
        log.debug(f"Finished on {x} {y}")
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

        val = 0.0
        dummy = 0.0
        valdx = 0.0
        valdy = 0.0
        valdxdy = 0.0
        valdxdx = 0.0
        valdydy = 0.0
        self.c_bicubic.getAllValues(x, y, val, dummy, dummy, valdxdy)
        self.c_bicubic.getSecondDerivatives(x, y, valdxdx, valdydy)

        return valdxdx * valdydy - valdxdy * valdxdy

    def getPoints(self):
        return self._store