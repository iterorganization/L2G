"""This module contains data container classes used to collect data from FLT
runs.
"""
import numpy as np
import numpy.typing as npt
import logging
import typing

log = logging.getLogger(__name__)

if typing.TYPE_CHECKING:
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
    from vtkmodules.vtkCommonDataModel import vtkPolyData


class L2GResults:
    """This class holds FLT relevant data.

    Size of attributes:
     * If 1D - the size corresponds to the number of target cells
     * if 2D - the size corresponds to the number of target cells * 3 (vector)

    Attributes:
        baryCent: 1D array of barycenters of triangular target input mesh
        normals: 2D Array of vectors of normals on baryCent
        flux: 1D Array of Poloidal magnetic flux in Webb/rad on baryCent
        BVec: 2D array of magnetic field vector in Cartesian coordinate
              system in T on baryCent
        BVecCyln: 2D array of the magnetic field vector in Cylndrical
                  coordinate system in T on baryCent. In this case we only have
                  poloidal and toroidal component.
        Bdot: 1D array Dot product between normals and BVec on baryCent
        angle: 1D array Incident angle between normals and BVec on baryCent
        mask: 1D array Intersection mask (1 - shadowed, 0 - wetted)
        direction: 1D array FL direction from a triangle
        conlen: 1D array Calculated connection length on a triangle in meters
        drsep: 1D array of radial distance along the midplane on baryCent in meters
        q: 1D array of incident heat load values on baryCent in W
        qpar: 1D array of parallel heat laod values on baryCent in W
        geom_hit_ids: 1D array of integer IDs of geometries loaded into TLAS
        vertices: 1D array of mesh vertices (size is equal to 3*N, where N is the
                number of vertices)
        triangles: 1D array of mesh triangles (size is equal to 3*M, where M is
                then number of triangles)

    """

    __slots__: list[str] = [
        "normals",  # Provided by FieldLineTracer
        "flux",     # Provided by FieldLineTracer
        "BVec",     # Provided by FieldLineTracer
        "BVecCyln", # Provided by FieldLineTracer
        "Bdot",     # Provided by FieldLineTracer
        "angle",    # Provided by FieldLineTracer
        "mask",     # Provided by FieldLineTracer
        "direction",# Provided by FieldLineTracer
        "conlen",   # Provided by External C++
        "drsep",    # Provided by FieldLineTracer
        "drsep2",   # Provided by FieldLineTracer
        "geom_hit_ids", # Provided by External C++
        "prim_hit_ids", # Provided by External C++
        "empty",
        "recalculate", # Used by processMagneticData in FieldlineTracer.
        "triangles",
        "vertices"
        ]

    arrays_to_dump: list[str] = [
        "normals",
        "flux",
        "BVec",
        "BVecCyln",
        "Bdot",
        "angle",
        "mask",
        "direction",
        "conlen",
        "drsep",
        "drsep2",
        "geom_hit_ids",
        "prim_hit_ids",
        ]

    def __init__(self):
        self.reset()

    def reset(self, N: int=0):
        if N == 0:
            self.normals:       npt.NDArray[np.float64] = np.array([])
            self.flux:          npt.NDArray[np.float64] = np.array([])
            self.BVec:          npt.NDArray[np.float64] = np.array([])
            self.BVecCyln:      npt.NDArray[np.float64] = np.array([])
            self.Bdot:          npt.NDArray[np.float64] = np.array([])
            self.angle:         npt.NDArray[np.float64] = np.array([])
            self.mask:          npt.NDArray[np.float64] = np.array([])
            self.direction:     npt.NDArray[np.int32] = np.array([])
            self.conlen:        npt.NDArray[np.float64] = np.array([])
            self.drsep:         npt.NDArray[np.float64] = np.array([])
            self.drsep2:        npt.NDArray[np.float64] = np.array([])
            self.geom_hit_ids:  npt.NDArray[np.int32] = np.array([])
            self.prim_hit_ids:  npt.NDArray[np.int32] = np.array([])

        else:
            # First we check if we already have arrays of the given size

            if self.flux is not None and self.flux.size == N:
                # Do nothing. Maybe we should at least zero the arrays.
                pass
            else:
                self.normals      = np.empty(3*N, dtype=np.float64)
                self.flux         = np.empty(N, dtype=np.float64)
                self.BVec         = np.empty(3*N, dtype=np.float64)
                self.BVecCyln     = np.empty(2*N, dtype=np.float64)
                self.Bdot         = np.empty(N, dtype=np.float64)
                self.angle        = np.empty(N, dtype=np.float64)
                self.mask         = np.empty(N, dtype=np.float64)
                self.direction    = np.empty(N, dtype=np.int32)
                self.conlen       = np.empty(N, dtype=np.float64)
                self.drsep        = np.empty(N, dtype=np.float64)
                self.drsep2       = np.empty(N, dtype=np.float64)
                self.geom_hit_ids = np.empty(N, dtype=np.int32)
                self.prim_hit_ids = np.empty(N, dtype=np.int32)

        self.empty: bool = True
        self.recalculate: bool = True

    def generate_vtk_object(self) -> 'vtkUnstructuredGrid':
        """From the data it generates a VTK object.
        """
        if self.triangles is None or self.vertices is None:
            log.error("No data to dump")
            return

        try:
            from vtkmodules.vtkCommonCore import vtkPoints
            from vtkmodules.vtkCommonDataModel import (vtkUnstructuredGrid,
                                                       vtkCellArray,
                                                       vtkTriangle)
            from vtkmodules.util.vtkConstants import VTK_TRIANGLE
            from vtkmodules.util import numpy_support
        except ImportError as e:
            log.error("VTK not available")
            return None

        data = vtkUnstructuredGrid()
        points = vtkPoints()
        cells = vtkCellArray()
        N_vertices = len(self.vertices) // 3
        for i in range(N_vertices):
            points.InsertPoint(i, self.vertices[i*3:i*3+3])
        N_cells = len(self.triangles) // 3
        for i in range(N_cells):
            triangle = vtkTriangle()
            pointIds = triangle.GetPointIds()

            pointIds.SetId(0, self.triangles[3*i])
            pointIds.SetId(1, self.triangles[3*i+1])
            pointIds.SetId(2, self.triangles[3*i+2])
            cells.InsertNextCell(triangle)
        data.SetPoints(points)
        data.SetCells(VTK_TRIANGLE, cells)

        cellData = data.GetCellData()
        for key in self.arrays_to_dump:
            array = getattr(self, key)
            if array is None:
                continue

            if key == "angle":
                # Transform the angle
                array = np.rad2deg(array)
                array = np.where(array > 90.0, array - 90.0, 90.0 - array)

            if key in ["drsep", "drsep2", "conlen"]:
                # Scale to mm from m
                array = np.array(array, copy=True) * 1000

            if len(array) == N_cells * 3:
                # Vector field
                tmp_buff = np.asarray(array)
                scalarField = numpy_support.numpy_to_vtk(tmp_buff.reshape(N_cells,

                                                         3))
            else:
                scalarField = numpy_support.numpy_to_vtk(array)
            scalarField.SetName(key)
            cellData.AddArray(scalarField)
        return data

class L2GFLs:
    """Class that holds calculated FLs on given target triangles.
    """

    __slots__: list[str] = ["points", "target_to_m"]

    def __init__(self):
        # List of Ids
        self.points: np.ndarray | None = None
        self.target_to_m: float = 1.0 # What inverse scale to use to
                                                # transform to the same
                                                # unit dimension as used
                                                # target.

    def generate_vtk_object(self) -> 'vtkPolyData':
        if self.points is None:
            raise Exception("No data to dump")

        try:
            from vtkmodules.vtkCommonCore import vtkPoints
            from vtkmodules.vtkCommonDataModel import (vtkCellArray,
                                                       vtkPolyLine,
                                                       vtkPolyData)
            from vtkmodules.util import numpy_support
        except ImportError as _:
            raise Exception("VTK not available")

        N_points = np.sum([len(_) for _ in self.points]) // 3

        data = vtkPolyData()
        points = vtkPoints()
        cells = vtkCellArray()
        points.SetNumberOfPoints(N_points)
        id_offset = 0

        def getLineLength(x, y, z):
            dx = np.diff(x)
            dy = np.diff(y)
            dz = np.diff(z)
            length = np.sqrt(dx * dx + dy * dy + dz * dz).sum()
            return length

        # Add colors to the lines
        lengths = []

        for line in self.points:
            r = np.asarray(line[::3]) / self.target_to_m
            z = np.asarray(line[1::3]) / self.target_to_m
            phi = np.asarray(line[2::3])
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            nfl = len(r)
            lengths.append(getLineLength(x, y, z))
            polyline = vtkPolyLine()
            polylineIds = polyline.GetPointIds()
            polylineIds.SetNumberOfIds(nfl)
            for i in range(nfl):
                points.SetPoint(id_offset + i, [x[i], y[i], z[i]])
                polylineIds.SetId(i, id_offset + i)
            cells.InsertNextCell(polyline)
            id_offset += nfl

        data.SetPoints(points)
        data.SetLines(cells)

        cellData = data.GetCellData()
        conlenScalar = numpy_support.numpy_to_vtk(lengths)
        conlenScalar.SetName('conlen')
        cellData.AddArray(conlenScalar)
        return data

class L2GResultsHLM:
    """Class for holding HLM results

    Additional arrays holds additional arrays, such as in case of Flat-Top and
    Ramp-DOwn where we also have conservative values.


    Attributes:
        q_inc: Perpendicular shadowed heat load
        q_par: Parallel unshadowed heat load
        flux_expansion: Flux expansion
        additional_arrays: Storage for additional arrays
    """
    __slots__: list[str] = ["q_inc", "q_par", "flux_expansion",
                            "additional_arrays", "empty"]

    def __init__(self):
        self.reset()

    def reset(self):
        self.q_inc: np.ndarray = np.array([])
        self.q_par: np.ndarray = np.array([])
        self.flux_expansion: np.ndarray = np.array([])
        self.additional_arrays: list[np.ndarray] = []
        self.empty: bool = True
