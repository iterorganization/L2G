from vtkmodules.vtkIOXML import (vtkXMLUnstructuredGridReader,
                                 vtkXMLUnstructuredGridWriter)
from vtkmodules.vtkIOLegacy import (vtkUnstructuredGridReader,
                                    vtkUnstructuredGridWriter)
from vtkmodules.util import numpy_support
from vtkmodules.util.vtkConstants import VTK_TRIANGLE
from vtkmodules.vtkFiltersExtraction import vtkExtractCellsByType
from vtkmodules.vtkCommonDataModel import (vtkUnstructuredGrid, vtkCellArray,
                                           vtkTriangle)
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter

import numpy as np
import os
import logging
log = logging.getLogger(__name__)

__readers__ = {}
__writers__ = {}
__vtk_obj__ = {}

class WrongVTKFileExtension(Exception):
    pass

class CantReadVTKFile(Exception):
    pass

def getReader(file_path: str) -> vtkXMLUnstructuredGridReader | vtkUnstructuredGridReader:
    global __readers__
    if file_path in __readers__:
        return __readers__[file_path]

    file_ext = os.path.splitext(file_path)[-1]
    if file_ext.startswith("."):
        file_ext = file_ext[1:]
    reader: vtkUnstructuredGridReader | vtkXMLUnstructuredGridReader
    if file_ext == "vtu":
        reader = vtkXMLUnstructuredGridReader()
    elif file_ext == "vtk":
        reader = vtkUnstructuredGridReader()
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.ReadAllNormalsOn()
        reader.ReadAllTensorsOn()
        reader.ReadAllTCoordsOn()
        reader.ReadAllFieldsOn()
    else:
        raise WrongVTKFileExtension
    reader.SetFileName(file_path)
    __readers__[file_path] = reader
    return reader

def getWriter(file_path: str) -> vtkXMLUnstructuredGridWriter | vtkUnstructuredGridWriter:
    """Support only VTUs!!!
    """
    global __writers__
    if file_path in __writers__:
        return __writers__[file_path]
    file_ext = os.path.splitext(file_path)[-1]
    if file_ext.startswith("."):
        file_ext = file_ext[1:]

    if file_ext == "vtu":
        writer = vtkXMLUnstructuredGridWriter()
    elif file_ext == "vtk":
        writer = vtkUnstructuredGridWriter()
    else:
        raise WrongVTKFileExtension

    __writers__[file_path] = writer
    return writer

def openFile(*args, **kwargs) -> None:
    """VTK does not open file. At least the VTK legacy and VTU files are not
    "open".
    """
    global __readers__, __writers__
    pass

def getMeshData(file_path: str) -> tuple[np.ndarray, np.ndarray]:
    """Get points and triangles from vtkUnstructuredGrid.
    """
    reader = getReader(file_path)
    if reader is None:
        raise CantReadVTKFile(f"Could not read VTK file {file_path}")

    reader.Update()
    obj = reader.GetOutput()
    vertex = numpy_support.vtk_to_numpy(obj.GetPoints().GetData())

    # Filter the celly by type
    filter_obj = vtkExtractCellsByType()
    filter_obj.SetInputData(obj)
    filter_obj.AddCellType(VTK_TRIANGLE)
    filter_obj.Update()

    out_obj = filter_obj.GetOutput()

    cells = out_obj.GetCells()
    if cells is None:
        return vertex, np.array([], dtype=np.uint32)

    c_data = numpy_support.vtk_to_numpy(cells.GetData())

    # The cells data holds information about one cell in the following way:
    # Cell_id id1 id2 id3. Where cell_id is set to the triangles id.
    # Therefore to get the ids we need to copy the data without the Cell
    # IDs.
    triangles = np.zeros((c_data.shape[0] // 4, 3), dtype=np.uint32)
    triangles[:, 0] = c_data[1::4]
    triangles[:, 1] = c_data[2::4]
    triangles[:, 2] = c_data[3::4]
    print(triangles)
    return vertex, triangles

def generateVtkObject(vertices, triangles):
    """Populate a vtkUnstructuredGrid object with the points and triangles
    """
    if triangles is None or vertices is None:
        return

    vtk_obj = vtkUnstructuredGrid()
    points = vtkPoints()
    cells = vtkCellArray()
    n_vertices = len(vertices)
    for i in range(n_vertices):
        points.InsertPoint(i, vertices[i])
    n_cells = len(triangles)
    for i in range(n_cells):
        triangle = vtkTriangle()
        pointIds = triangle.GetPointIds()
        pointIds.SetId(0, triangles[i][0])
        pointIds.SetId(1, triangles[i][1])
        pointIds.SetId(2, triangles[i][2])
        cells.InsertNextCell(triangle)
    vtk_obj.SetPoints(points)
    vtk_obj.SetCells(VTK_TRIANGLE, cells)
    return vtk_obj

def getVtkObj(file_path: str, vertices: np.ndarray, triangles: np.ndarray):
    global __vtk_obj__
    """Holding a reference VTK object on a file_path.
    """


    if not file_path in __vtk_obj__:
        vtk_obj = generateVtkObject(vertices, triangles)
        __vtk_obj__[file_path] = vtk_obj

    return __vtk_obj__[file_path]

def writeMesh(file_path, vertices, triangles, mesh_name: str='mesh'):
    """Depending on the file_path extension we create a writer.
    """

    writer = getWriter(file_path)
    if writer is None:
        log.error(f"Could not get a writer for the {file_path}")
        return
    vtk_obj = getVtkObj(file_path, vertices ,triangles)

    writer.SetFileName(file_path)
    writer.SetInputData(vtk_obj)
    writer.Update()

def getNumberOfTimeSteps(file_path: str, *args) -> int:
    """Return the number of time steps from inside the VTK file. VTK legacy
    file can only have data of one time step.

    Returns:
        number_of_time_steps (int): 1 if Legacy VTK file, else however number
            of time steps are in a VTU file.
    """
    reader = getReader(file_path)
    if file_path.endswith("vtk") or isinstance(reader, vtkUnstructuredGridReader):
        return 1

    reader.Update()
    return reader.GetNumberOfTimeSteps()

def addArrayToVtkObj(vtk_obj: vtkUnstructuredGrid, array_name: str,
                     array: np.ndarray):
    """To an input vtk_obj add an array array with name array_name.

    Arguments:
        vtk_obj (vtkUnstructuredGrid): VTK obj
        array_name (str): Name of the array
        array (np.ndarray): Numpy 1D array.
    """
    cell_data = vtk_obj.GetCellData()

    vtk_array = numpy_support.numpy_to_vtk(array)
    vtk_array.SetName(array_name)
    cell_data.AddArray(vtk_array)

def copyVtkObj(orig_obj: vtkUnstructuredGrid) -> vtkUnstructuredGrid:
    """Creates a deep copy of a provided VTK object.

    Arguments:
        orig_obj (vtkUnstructuredGrid): Input VTK object to copy

    Returns:
        obj (vtkUnstructuredGrid): Deeply copied VTK object.
    """
    obj = vtkUnstructuredGrid()
    obj.DeepCopy(orig_obj)
    return obj

def checkIfFieldExists(file_path: str, field_name: str) -> bool:
    # Do nothing for now
    return True

def getAllFieldIterations(file_path: str, field_name: str) -> list:
    return []

def doesItContainGroup(*args, **kwargs) -> bool:
    return False

def getGroupArr(*args, **kwargs) -> np.ndarray:
    return np.array([])

def getAllGroups(*args, **kwargs) -> list:
    return []

def getAllFieldNames(file_path: str) -> list[str]:
    """Return the name of all fields.

    Argument:
        file_path (str): Path to the VTK file.

    Returns:
        list (list): List of all field names.
    """

    reader = getReader(file_path)
    reader.Update()

    obj = reader.GetOutput()
    cell_data = obj.GetCellData()
    n_arrays = cell_data.GetNumberOfArrays()
    return [cell_data.GetArrayName(i) for i in range(n_arrays)]

def getField(file_path: str, field_name: str, index: int) -> np.ndarray:
    """Get the field, or array in VTK terms, from a VTK file.

    Arguments:
        file_path (str): Path to VTK file
        field_name (str): Name of the VTK array
        index (int): The index of the time step.

    Returns:
        array (np.ndarray): Numpy 1D array.
    """
    reader = getReader(file_path)
    if not isinstance(reader, vtkUnstructuredGridReader):
        reader.SetTimeStep(index)
    reader.Update()

    obj = reader.GetOutput()
    cell_data = obj.GetCellData()
    n_arrays = cell_data.GetNumberOfArrays()
    for i in range(n_arrays):
        name = cell_data.GetArrayName(i)
        if name == field_name:
            array = numpy_support.vtk_to_numpy(cell_data.GetArray(i))
            return array
    return np.array([])

def getMeasurements(file_path: str) -> np.ndarray:
    """Get cell areas (of triangular cells) from a VTK file.

    Arguments:
        file_path (str): Path to VTK file.

    Returns:
        array (np.ndarray): 1D array containing cell areas.
    """
    reader = getReader(file_path)
    if reader is None:
        log.error(f'Could not read this VTK file {file_path}')
        return np.array([])
    reader.Update()

    cell_size_filter = vtkCellSizeFilter()
    cell_size_filter.SetInputData(reader.GetOutput())
    cell_size_filter.SetComputeArea(True)
    cell_size_filter.Update()

    output = cell_size_filter.GetOutput()
    cell_data = output.GetCellData()

    array = numpy_support.vtk_to_numpy(cell_data.GetArray('Area'))
    return array
