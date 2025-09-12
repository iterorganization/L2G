import l2g.external.bicubic
import numpy as np

class PyRKF45FLT:
    """Class that wraps a C++ implementation of the RKF4(5) method for solving
    the fieldline tracing equations in an axisymmetrical Cylindrical coordinate
    system.

    .. note::

       This is not an implementation that can solve a general system of PDEs.
       This class is meant for easily plotting specific fieldlines.

    In other words it is a class that follows contours of a 2D function,
    defined in a (R, Z) Cylindrical space. Based on the

    .. code_block:: python

       from l2g.external.bicubic import PyBicubic
       from l2g.external.rkf45 import PyRKF45FLT


       int_obj = PyBicubic(x, y, f)

       rkf45_obj = PyRKF45FLT()

       # Set the contour data
       rkf45_obj.setInterpolator(int_obj)
       # Essentially set the constant value F = B_t * R, from which B_t is
       # evaluated as F / R.
       rkf45_obj.set_vacuum_fpol(vacuum_fpol)


       # Start following
       r: float
       z: float
       th: float
       time_step: float
       nr: float
       nz: float
       nth: float

       while 1:
           nr, nz, nth = rkf45_obj.run_step(r, z, th, time_step)
    """
    def __init__(self) -> None: ...
    def setInterpolator(self, obj: l2g.external.bicubic.PyBicubic) -> None:
        """Sets the interpolator object that contains the data for following
        contours.

        Arguments:
            obj (PyBicubic): Bicubic interpolator used for obtaining the
                gradient values of the interpolated data.
        """
    def set_r_move(self, r_move: float) -> None:
        """Sets the radial shift of the data.

        Arguments:
            r_move (float): Radial shift.
        """
    def set_z_move(self, z_move: float) -> None:
        """Sets the vertical shift of the data.

        Arguments:
            r_move (float): Vertical shift.
        """
    def set_vacuum_fpol(self, vacuum_fpol: float) -> None:
        """Sets the constant F = B_t R, used for calculating the toroidal
        component of the magnetic field as B_t = F / R.

        Arguments:
            vacuum_fpol (float): Poloidal current value.
        """
    def run_step(self, r: float, z: float, th: float, time_step: float) -> tuple[float, float, float]:
        """Perform one step of RKF4(5) from point (r, z, th) in the Cylindrical
        coordinate system.

        Arguments:
            r (float): Radial position
            z (float): Vertical position
            th (float): Toroidal angle
            time_step (float): Size of time_step

        Returns:
            out (tuple[float, float, float]): nr, nz, nth
        """
    def run_n_steps(self, r: float, z: float, th: float, time_step: float, n: int) -> tuple[float, float, float]:
        """Perform n steps of RKF4(5) from point (r, z, th) in the Cylindrical
        coordinate system.

        Arguments:
            r (float): Radial position
            z (float): Vertical position
            th (float): Toroidal angle
            time_step (float): Size of time_step
            n (int): Number of steps of size time_step to perform

        Returns:
            out (tuple[float, float, float]): nr, nz, nth
        """
