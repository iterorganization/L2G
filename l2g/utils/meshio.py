import numpy as np
import os
from typing import Tuple

import medcoupling as mc

import logging
log = logging.getLogger(__name__)

"""Contains functions for manipulating MED fields in terms of transforming
from MED Field to numpy array, vice versa and to write MED Field to existing
mesh.
"""

class FieldComponentException(Exception):
    """Exception that is raised, when trying to convert a numpy array to MED
    field. If we provide an array which is N-dimensional and provide a
    component info list which is of size M and if N != M, then the exception is
    raised, as the number of components (or labels for dimensions) has to be
    the same as the number of dimensions.
    """
    def __init__(self,
                 msg="Mismatch between shape of array and size of component"):
        super(FieldComponentException, self).__init__(msg)

def fieldToNumpy(file: str, field_name: str, iteration: int = 0, order: int = -1
    ) -> np.ndarray:
    """This function is used to read a field from a MED file.

    Returns:
        array (np.ndarray): Numpy array of the requested field

    Arguments:
        file (str): Path to the MED file
        field_name (str): Name of the requested field.
        iteration (int): Slice iteration. Default 0 or first
        order (int): Secondary slice iteration. Not used currently,
            so default is always -1

    """
    field = mc.ReadField(file, field_name, iteration, order)
    array = field.getArray()
    return array.toNumPyArray()

def numpyArrayToField(arr: np.ndarray, field_name: str,
        mesh: mc.MEDCouplingUMesh, associated_time: int = 0, iteration: int = 0,
        info_on_component: list=[]) -> mc.MEDCouplingFieldDouble:
    """Converts a numpy array to field.

    Returns:
        field (mc.MEDCouplingFieldDouble): Prepared MED Field to be
            written/added to a MED file.

    Arguments:
        arr (np.ndarray): Numpy array to convert to MED Field
        field_name (str): Name of field
        mesh (mc.MEDCouplingUMesh): MED mesh
        associated_time (int): Associated time for the field. In seconds
        iteration (int): Iteration of the Field, not to be confused with
            associated time. Iteration merely specifies the index of the field
            in a MED file (0 start index).
        info_on_component (list): When we have a component field, e.g., a N
            dimensional vector field, then this is how we tell this function to
            properly re-shape the numpy array so that when we view the data in
            ParaView, the data is properly labeled and displayed.

    """
    # Get shape for Array component! Meaning if it is a vector then the array
    # needs to know that


    # Find out how many values per cells there are
    if info_on_component:
        # Vector field. Check if the length of the component info is the same
        # as the ratio between the number of values in a numpy array and.
        num_of_components = len(info_on_component)
        N_cells = mesh.getNumberOfCells()
        ratio = int(len(arr) / N_cells)
        if ratio == num_of_components:
            arr = arr.reshape((N_cells, ratio))
        else:
            # Raise exception, otherwise medcoupling will raise it due to
            # mismatch in numbers of components and the size of the field. Or
            # not, maybe it will trigger just an UB.
            msg = f"Shape of array {arr.shape} does not match the number of "
            msg += f"components {info_on_component}"
            raise FieldComponentException(msg)

    mcArray = mc.DataArrayDouble(np.asarray(arr, dtype=np.float))
    if info_on_component:
        mcArray.setInfoOnComponents(info_on_component)

    field = mc.MEDCouplingFieldDouble.New(mc.ON_CELLS, mc.ONE_TIME)
    field.setMesh(mesh)
    field.setTime(associated_time, iteration, -1)
    field.setTimeUnit("s")
    field.setName(field_name)
    field.setArray(mcArray)
    return field

def writeFieldToAlreadyExistingMesh(field: mc.MEDCouplingFieldDouble,
    med_file: str) -> None:
    """Writes a MED Field to a MED file.

    Arguments:
        field (mc.MEDCouplingFieldDouble): Field to be written to a MED file
        med_file (str): Path to the MED file
    """
    fMEDFile=mc.MEDFileField1TS.New()
    fMEDFile.setFieldNoProfileSBT(field)
    fMEDFile.write(med_file,0)
    # mc.WriteFieldUsingAlreadyWrittenMesh(med_file, field)

def readVtkMesh(filePath: str) -> Tuple[np.array, np.ndarray]:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(filePath)
    # Irrelevant as we need only nodes and cells.
    #reader.ReadAllVectorsOn()
    #reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()
    n_points = data.GetNumberOfPoints()
    nodes = vtk_to_numpy(data.GetPoints().GetData())
    nodes = nodes.reshape(n_points * 3)

    n_cells = data.GetNumberOfCells()
    cells = np.empty(n_cells * 3, dtype=np.uint32)
    _c = 0

    # Process cells... one by one, as there is no GetCellsByType function...
    for i in range(n_cells):
        cell = data.GetCell(i)
        if isinstance(cell, vtk.vtkTriangle):
            cellId = cell.GetPointIds()
            cells[3 * _c] = cellId.GetId(0)
            cells[3 * _c + 1] = cellId.GetId(1)
            cells[3 * _c + 2] = cellId.GetId(2)

            _c += 1

    # Remove unpopulated values
    if _c != n_cells:
        cells = cells[:_c]

    return nodes, cells

def readMedMesh(med_path: str) -> Tuple[np.array, np.ndarray]:
    """Gets triangular mesh from a MED file. The output is a tuple that
    contains 1D array of vertices (point data flatten to 1D) and 1D array of
    triangles (triangle ID data flatten to 1D).

    Returns:
        tuple (tuple): A tuple consisting of 1D array of vertices and 1D array
            of triangles.

    Arguments:
        med_path (str): Path  to MED file
    """
    mesh = mc.ReadMeshFromFile(med_path)

    # Vertexes
    vertex = mesh.getCoords().toNumPyArray()
    vertex = np.array(vertex, dtype=np.float32)

    # Triangles

    Ntriangles = mesh.getNumberOfCellsWithType(mc.NORM_TRI3)

    triangles = np.empty(Ntriangles * 3, np.uint32)

    Ncells = mesh.getNumberOfCells()

    c = 0
    for i in range(Ncells):
        if mesh.getTypeOfCell(i) != mc.NORM_TRI3:
            continue

        triangles[c:c+3] = mesh.getNodeIdsOfCell(i)
        c += 3

    return vertex.reshape(vertex.shape[0] * 3), triangles

def readMesh(filePath: str) -> Tuple[np.ndarray, np.ndarray]:
    if not os.access(filePath, os.R_OK):
        log.error(f"Cannot read file {filePath}!")
        raise IOError

    ext = filePath.split('.')[-1].lower()

    if ext == 'med':
        log.debug(f"Reading MED {filePath}.")
        return readMedMesh(filePath)
    elif ext == 'vtk':
        log.debug(f"Reading VTK {filePath}.")
        return readVtkMesh(filePath)
    else:
        log.error(f"Unknown extension: {ext}!")
        return None, None
