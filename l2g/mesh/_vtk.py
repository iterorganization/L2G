from typing import List, Tuple
import vtk
from vtk.util import numpy_support
import numpy as np
import os
import logging
log = logging.getLogger(__name__)

__readers__ = {}
__writers__ = {}
__vtk_obj__ = {}

def getReader(file_path: str):
    global __readers__
    if file_path in __readers__:
        return __readers__[file_path]

    file_ext = os.path.splitext(file_path)[-1]
    if file_ext.startswith("."):
        file_ext = file_ext[1:]
    reader = None
    if file_ext == "vtu":
        from vtk import vtkXMLUnstructuredGridReader
        reader = vtkXMLUnstructuredGridReader()
    elif file_ext == "vtk":
        from vtk import vtkUnstructuredGridReader
        reader = vtkUnstructuredGridReader()
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.ReadAllNormalsOn()
        reader.ReadAllTensorsOn()
        reader.ReadAllTCoordsOn()
        reader.ReadAllFieldsOn()
    reader.SetFileName(file_path)
    __readers__[file_path] = reader
    return reader

def getWriter(file_path: str):
    """Support only VTUs!!!
    """
    global __writers__
    if file_path in __writers__:
        return __writers__[file_path]
    file_ext = os.path.splitext(file_path)[-1]
    if file_ext.startswith("."):
        file_ext = file_ext[1:]
    writer = None
    if file_ext == "vtu":
        from vtk import vtkXMLUnstructuredGridWriter
        writer = vtkXMLUnstructuredGridWriter()
    elif file_ext == "vtk":
        from vtk import vtkUnstructuredGridWriter
        writer = vtkUnstructuredGridWriter()
        writer.ReadAllScalarsOn()
        writer.ReadAllVectorsOn()
        writer.ReadAllNormalsOn()
        writer.ReadAllTensorsOn()
        writer.ReadAllTCoordsOn()
        writer.ReadAllFieldsOn()
    __writers__[file_path] = writer
    return writer

def openFile(file_path, **kwargs) -> None:
    global __readers__, __writers__

def getMeshData(file_path: str) -> Tuple[np.ndarray]:
    """Get points and triangles from vtkUnstructuredGrid.
    """
    reader = getReader(file_path)
    if reader is None:
        log.error(f"Could not read this VTK file {file_path}")
        return
    reader.Update()
    obj = reader.GetOutput()
    vertex = numpy_support.vtk_to_numpy(obj.GetPoints().GetData())

    # Filter the celly by type
    filter_obj = vtk.vtkExtractCellsByType()
    filter_obj.SetInputData(obj)
    filter_obj.AddCellType(vtk.VTK_TRIANGLE)
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

    vtk_obj = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    n_vertices = len(vertices)
    for i in range(n_vertices):
        points.InsertPoint(i, vertices[i])
    n_cells = len(triangles)
    for i in range(n_cells):
        triangle = vtk.vtkTriangle()
        pointIds = triangle.GetPointIds()
        pointIds.SetId(0, triangles[i][0])
        pointIds.SetId(1, triangles[i][1])
        pointIds.SetId(2, triangles[i][2])
        cells.InsertNextCell(triangle)
    vtk_obj.SetPoints(points)
    vtk_obj.SetCells(vtk.VTK_TRIANGLE, cells)
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

def getNumberOfTimeSteps(file_path: str):
    if file_path.endswith("vtk"):
        return 1

    reader = getReader(file_path)
    reader.Update()
    return reader.GetNumberOfTimeSteps()

def addArrayToVtkObj(vtk_obj: vtk.vtkUnstructuredGrid, array_name: str,
                     array: np.ndarray):
    cell_data = vtk_obj.GetCellData()

    vtk_array = numpy_support.numpy_to_vtk(array)
    vtk_array.SetName(array_name)
    cell_data.AddArray(vtk_array)

def copyVtkObj(orig_obj: vtk.vtkUnstructuredGrid):
    obj = vtk.vtkUnstructuredGrid()
    obj.DeepCopy(orig_obj)
    return obj

def checkIfFieldExists(file_path: str, field_name: str) -> bool:
    # Do nothing for now
    return True

def getAllFieldIterations(file_path: str, field_name: str) -> List:
    return []

def getAllFieldNames(file_path: str):
    return []

def getField(file_path: str, field_name: str, index: int) -> np.ndarray:
    reader = getReader(file_path)
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
