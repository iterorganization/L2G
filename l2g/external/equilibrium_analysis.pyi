import numpy as np
import l2g.equil

class EQA:
    equilibrium: l2g.equil.Equilibrium
    psi_grad_sign: float
    evaluated: bool
    lcfs_points: list

    def __init__(self, equilibrium: l2g.equil.Equilibrium | None = None) -> None: ...
    def resetValues(self) -> None: ...
    def setEquilibrium(self, obj: l2g.equil.Equilibrium) -> None:
        """Sets the equilibrium objects and propagates the
        flux data to the interpolation and solver objects.
        """
    def setDisplacement(self, r_displ: float, z_displ: float) -> None: ...
    def evaluate(self) -> None:
        """Evaluates the equilibrium, determining the type (div or lim) and
        obtaining characteristic data.
        """
    def getType(self) -> str:
        """Plasma type.
        Returns:
            plasma_type (str): Either 'div' or 'lim'
        """
    def getBoundaryFluxValue(self) -> float:
        """Poloidal flux value of the boundary magnetic surface (either a
        closed surface in a limiter configuration or a separatrix during
        diverted operation). In Webb/rad.

        Retruns:
            psi (float): Flux value. In Webb/rad.
        """
    def getSecondaryXFluxValue(self) -> float:
        """Poloidal flux value of the secondary separatrix, during a double-X
        diverted operation. In Webb/rad.

        Returns:
            psi (float): Flux value. In Webb/rad.
        """
    def getContactPoint(self) -> list[float]:
        """Get the contact point (R, Z) on the wall of a limiter operation.

        Returns:
            contact_point (list): (R, Z) location.
        """
    def getLowerXPoint(self) -> list[float]:
        """Get the lower X point (R, Z) of a diverted operation.

        Returns:
            contact_point (list): (R, Z) location.
        """
    def getUpperXPoint(self) -> list[float]:
        """Get the upper X point (R, Z) of a diverted operation.

        Returns:
            contact_point (list): (R, Z) location.
        """
    def getMidplaneInfo(self, lcfs: float | None = None, which: str = "owl") -> tuple[float, float, float, float]:
        """Obtain the midplane information, mainly the point (Rb, Z) of where
        the midplane starts, B total [T] (magnitude of the total magnetic
        field at (Rb, Z)) and B pm (poloidal component of the magnetic field at
        (Rb, Z)).

        Arguments:
            lcfs (float | None): Magnetic surface magnetic poloidal flux value
                from which the midplane starts. Optional, by default the LCFS
                flux value is taken.
            which (str): Select which midplane, either inner 'iwl' or outer
                'owl'

        Returns:
            values (tuple[float, float, float, float]): Rb [m], Z[m],
                Btotal[T], Bpm[T]
        """
    def alingLcfsToPoints(self, point: list[float]) -> list[float]:
        """Finds the shortest displacement required to put the limiter LCFS
        to point p.

        Arguments:
            p (list): 1D array with 2 values.
        Returns:
            displ (np.ndarray): New evaluated displacement.
        """
    def distanceBetweenPsiOnMidplane(self, sep1: float, sep2: float) -> float:
        """Calculates the difference between the provided magnetic surface on
        the midplane.

        Returns:
            drsep (float): Distance between separatrixes, in mm
        """
    def getPsi(self, r: float, z: float) -> tuple[float, float, float]: ...
    def getPsiDxDy(self) -> np.ndarray: ...
