import l2g.external.bicubic

class PyBfgs2d:
    """Class that implements the BFGS (Broyden–Fletcher–Goldfarb–Shanno)
    algorithm for solving unconstrained nonlinear optimization problems.
    Particularly it uses the Bicubic interpolation algorithm of a 2D function
    and finds the minimum/maximum of the function in the definition area of
    the function.

    .. code-block::

       from l2g.external.bfgs_2d import PyBfgs2d
       from l2g.external.bicubic import PyBicubic

       int_obj = PyBicubic()

       # Load the input data
       ...

       opt_obj = PyBfgs2d()
       opt_obj.setInterpolator(int_obj)

       # Find the minimum
       guess_x: float
       guess_y: float
       max_iter: int


       out_of_bounds: bool
       x: float
       y: float
       out_of_bounds, x, y = obt_obj.findMinimum(guess_x, guess_y, max_iter)

       if not out_of_bounds:
           print("Success!")

    """
    def __init__(self) -> None: ...
    def setInterpolator(self, obj: l2g.external.bicubic.PyBicubic) -> None:
        """Sets the Bicubic interpolator.
        """
    def setBounds(self, x1: float, x2: float, y1: float, y2: float) -> None:
        """Set the bounds of the search area.

        Arguments:
            x1 (float): Left border of X interval.
            x2 (float): Right border of X interval.
            y1 (float): Left border of Y interval.
            y2 (float): Right border of Y interval.
        """
    def findMinimum(self, guess_x: float, guess_y: float, max_it: int) -> tuple[bool, float, float]:
        """This is a gradient-descent minimizer function that  tries to find a
        point in the function domain space where the gradient is close to
        zero.

        The main use of the function is finding the minimum of the function,
        hence the name remained, although it can get trapped and returns the
        position of a extrema point (saddle point).

        Arguments:
            guess_x (double): Starting X position
            gouss_y (double): Starting Y position
            max_it (int): Maximum number of iterations

        Returns:
            out_of_bounds (bool): True if the method went out of bounds. Is
                False when successful.
            x (double): X position of the solution.
            y (double): Y position of the solution.
        """
    def secondDerivativeTest(self, x: float, y: float) -> bool:
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

        Arguments:
            x (double): X position
            y (double): Y position
        Returns:
            D (double): Value of the second derivative test.

        """
    # For tracking progression during minimizing.
    def getPoints(self) -> list[list[int]]:
        """Get the points that tracks the progression of the algorithm from
        a starting point and to the solution.

        Returns:
            points (list): List of points of the form [[x1, y1], [x2, y2],...]
        """