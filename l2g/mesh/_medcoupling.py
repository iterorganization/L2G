from typing import Dict, Tuple
import medcoupling as mc
import numpy as np

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
        mesh: mc.MEDCouplingUMesh, associated_time: float = 0,
        iteration: int = 0, info_on_component: list=[]) -> mc.MEDCouplingFieldDouble:
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
    file_path: str) -> None:
    """Writes a MED Field to a MED file.

    Arguments:
        field (mc.MEDCouplingFieldDouble): Field to be written to a MED file
        file_path (str): Path to the MED file
    """
    fMEDFile=mc.MEDFileField1TS.New()
    fMEDFile.setFieldNoProfileSBT(field)
    fMEDFile.write(file_path,0)
    # mc.WriteFieldUsingAlreadyWrittenMesh(file_path, field)

__med__: Dict[str, mc.MEDCouplingUMesh] = {}
__file_med__: Dict[str, mc.MEDFileUMesh] = {}

def openFile(file_path: str, **kwargs) -> None:
    """Opens a file and creates MEDCoupling File and Object handles.

    Arguments:
        file_path (str): Path to MED file.
        kwargs (dict): Optional arguments passed. For instance a MED file can
            have multiple meshes, therefore we can also use the name of the
            mesh to obtain the desired mesh.
    """
    global __med__, __file_med__

    if 'name' in kwargs:
        __file_med__[file_path] = mc.MEDFileUMesh.New(file_path, kwargs['name'])
    else:
        __file_med__[file_path] = mc.MEDFileUMesh.New(file_path)
    __med__[file_path] = __file_med__[file_path].getMeshAtLevel(0)

def checkIfFieldExists(file_path: str, field_name: str) -> bool:
    """Checks if a field exists in the file.

    Arguments:
        file_path (str): Path to file
        field_name (str): Name of field

    Returns:
        ok (bool)
    """
    fields = mc.GetAllFieldNames(file_path)
    if field_name in fields:
        return True
    return False

def getMeshData(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Returns
    """
    global __med__, __file_med__
    if not file_path in __file_med__:
        openFile(file_path)
    med = __med__[file_path]
    vertices = np.array(med.getCoords().toNumPyArray(), np.float32)

    n_triangles = med.getNumberOfCellsWithType(mc.NORM_TRI3)
    triangles= np.empty((n_triangles, 3), np.uint32)
    c = 0
    for i in range(n_triangles):
        if med.getTypeOfCell(i) != mc.NORM_TRI3:
            continue
        triangles[c][:] = med.getNodeIdsOfCell(i)
        c += 1
    return vertices, triangles

def writeMesh(file_path, vertices, triangles, mesh_name: str="mesh"):
    """Saves an unstructured grid consisting of triangles to new file path.
    """

    # Create a new mesh file.

    med_umesh_file = mc.MEDFileUMesh.New()

    # Dimension level is 2
    med = mc.MEDCouplingUMesh.New(mesh_name, 2)
    # setCoords always expect DataArrayDouble
    med.setCoords(mc.DataArrayDouble(np.array(vertices, dtype=np.float64)))

    n_cells = triangles.shape[0]
    med.allocateCells(n_cells)

    for i in range(n_cells):
        med.insertNextCell(mc.NORM_TRI3, triangles[i].tolist())
    # med.checkConsistencyLight()

    # Set the mesh always at level 0.
    med_umesh_file.setMeshAtLevel(0, med)
    # Overwrite mesh at location.
    med_umesh_file.write(file_path, 0)

def getNumberOfTimeSteps(file_path: str, *args):
    if args:
        field_name = args[0]
    else:
        # Typical FLT field
        field_name = "conlen"

    iterations = mc.GetAllFieldIterations(file_path, field_name)
    return len(iterations)

def getField(file_path: str, field_name: str, index: int, order: int=-1):
    if checkIfFieldExists(file_path, field_name):
        return fieldToNumpy(file_path, field_name, index, order)
    else:
        return None

def writeField(file_path: str, array_name: str, array: np.ndarray, index: int,
               time: float, info_on_component: list = []):

    global __med__

    if file_path not in __med__:
        openFile(file_path)
    med = __med__[file_path]
    field = numpyArrayToField(arr=array, field_name=array_name, mesh=med, iteration=index,
            associated_time=time, info_on_component=info_on_component)
    writeFieldToAlreadyExistingMesh(field, file_path)