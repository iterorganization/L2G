import numpy as np

from typing import Optional

import logging
log = logging.getLogger(__name__)

class Equilibrium(object):
    """Python object class that holds Equilibrium data. It holds the minimum
    required data for FLT use. This also includes for applying hlm profiles.

    The minimum required data:

      - Wall silhouette. Used for determining from an equilibrium the type of
        the equilibrium
      - Magnetic axis location.
      - RZ grid and magnetic poloidal flux. For the FL tracing equations.
      - Ip Plasma current

    Attributes:

        wall_contour_r (arr): R (radial) points of the 2D wall silhouette. In
            meters
        wall_contour_z (arr): Z (vertical) points of the 2D wall silhouette. In
            meters

        grid_dim_r (arr): Number of points in the radial direction
        grid_dim_z (arr): Number of points in the vertical direction.
        grid_r (arr): 1D array of R (radial) points of the R-Z grid. In meters.
        grid_z (arr): 1D array of Z (vertical) points of the R-Z grid. In
            meters.
        psi (arr): 2D Array of magnetic poloidal flux values of size
            grid_dim_r x grid_dim_z. Row major ordered which means that rows,
            are concatenated in 1D array. In Webb/rad.
        psi_boundary (float): Value of psi at boundary.
        psi_axis (float): Value of psi at magnetic axis.
        psi_sign (int): Tells you the sign of the flux gradient when going
            away from plasma center to boundary.
        fpol_vacuum (float): FPol value (Bt R) in vacuum. In m T (meter Tesla).
        Ip (float): Plasma current.
        a (float): minor radius
        Psol (float): Power @ SOL
        Area (float): Plasma crossection area.
    """
    __slots__ = [
        "wall_contour_r",
        "wall_contour_z",
        "mag_axis_r",
        "mag_axis_z",
        "grid_dim_r",
        "grid_dim_z",
        "grid_r",
        "grid_z",
        "psi",
        "psi_boundary",
        "psi_axis",
        "Ip",
        "fpol",
        "fpol_flux",
        "fpol_vacuum",
        "psi_sign",
        "type",
        "a",
        "Psol",
        "Area"]

    def __init__(self):
        self.wall_contour_r: list = []
        self.wall_contour_z: list = []

        self.mag_axis_r: float = -1.0
        self.mag_axis_z: float = -1.0

        self.grid_dim_r: int = -1
        self.grid_dim_z: int = -1

        self.grid_r: np.ndarray = np.array([])
        self.grid_z: np.ndarray = np.array([])
        self.psi: np.ndarray = np.array([])
        self.psi_sign: Optional[int] = None # Sign of flux gradient from center to boundary.
        self.psi_boundary: Optional[float] = None
        self.psi_axis: Optional[float] = None

        self.Ip: Optional[float] = -1.0

        self.fpol: list = []
        self.fpol_flux: list = []
        self.fpol_vacuum: float= -1.0

        self.type: str = ''

        self.a: float = 0.0
        self.Psol: float = 0.0
        self.Area: float = 0.0

        pass

def correct_equilibrium_helicity(eq_obj: Equilibrium) -> None:
    """Corrects the equilibrium helicity. That is changes the sign of the
    poloidal flux values and corrects if necessary.

    In most convention we wish to have the gradient of the poloidal
    magnetic flux showing out.

    Relevant IDM document: GGKFNP
    """

    # if eq_obj.psi_axis is None:

    #     interp2d = SplineInterpolator()
    #     interp2d.interpolate_data(eq_obj.grid_r, eq_obj.grid_z, eq_obj.psi)

    #     # Get mag axis flux value
    #     mag_axis_psi = interp2d.ev(eq_obj.mag_axis_r, eq_obj.mag_axis_z)
    # else:
    #     mag_axis_psi = eq_obj.psi_axis

    # # Subtract from the boundary and see if it is falling or rising.
    # sign = np.sign(eq_obj.psi_boundary - mag_axis_psi)

    if eq_obj.psi_sign == -1:
        log.info("Incorrect Flux sign for ITER. Rectifying.")
        eq_obj.psi *= -1
        eq_obj.psi_boundary *= -1
        eq_obj.psi_axis *= -1
        eq_obj.psi_sign = 1
    else:
        log.info("Correct helicity for ITER.")

    # Check FPOL vacuum. Should show in the positive orientation or negative
    # direction.

    if eq_obj.fpol_vacuum > 0:
        log.info("Correcting FPOL (Bt R) to have negative sign (positive " +
                 "orientation) for ITER.")
        eq_obj.fpol_vacuum *= -1
        eq_obj.fpol *= -1
    else:
        log.info("Sign of FPOL (Bt R) correct.")

    return None
