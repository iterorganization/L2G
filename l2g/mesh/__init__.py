import numpy as np
import logging
import os
log = logging.getLogger(__name__)
from typing import Dict, Tuple, TYPE_CHECKING, Union

from l2g.comp import L2GResults, L2GResultsHLM, L2GFLs

def rotatePointsAroundAxis(points: np.ndarray, p1: np.ndarray, p2: np.ndarray,
    theta: float):
    """Rotate the argument points around the axis |p2 - p1| for the angle
    theta.

    Arguments:
        points (npt.ndarray): Array of points with the shape ((N_points, 3)).
        p1 (npt.ndarray): Start point of the axis
        p2 (npt.ndarray): End point of the axis
        theta (float): angle in degrees.

    Returns:
        rotated_points (np.ndarray)
    """

    # Copy the points and translate it to p1.
    out_points = np.zeros(points.shape)
    out_points[:] = points - p1

    # Get the normalized diff vector.
    diff_p = p2 - p1
    u = diff_p / np.linalg.norm(diff_p)

    # Construct the rotational matrix
    R = np.zeros((3, 3))

    cos_th = np.cos(np.deg2rad(theta))
    sin_th = np.sin(np.deg2rad(theta))
    moin_cos_th = 1 - cos_th
    ux = u[0]
    uy = u[1]
    uz = u[2]
    R[0, 0] = cos_th + ux * ux * moin_cos_th
    R[1, 0] = uy * ux * moin_cos_th + uz * sin_th
    R[2, 0] = uz * ux * moin_cos_th - uy * sin_th

    R[0, 1] = ux * uy * moin_cos_th - uz * sin_th
    R[1, 1] = cos_th + uy * uy * moin_cos_th
    R[2, 1] = uz * uy * moin_cos_th + ux * sin_th

    R[0, 2] = ux * uz * moin_cos_th + uy * sin_th
    R[1, 2] = uy * uz * moin_cos_th - ux * sin_th
    R[2, 2] = cos_th + uz * uz * moin_cos_th

    # Rotate the points
    # out_points = R@out_points.T
    out = (R@out_points.T).T

    # We need to fix the stride. Best way is to copy all values...
    # No idea otherwise how to fix the array stride.
    out_points[:] = out
    # out_points = out_points.T

    # Re-translate it to correct place.
    # out_points += p1
    return np.array(out_points + p1, dtype=np.float32)

class Mesh():
    """This is a convenient class for reading 2D surface meshes made of
    triangles. Additionally it writes the mesh in different profiles and writes
    fields to those files. It supports HDF5, MEDCoupling and VTK file formats.

    The idea is to have a way to write results in a file for storage.

    For visualization you can open files saved in the VTK file format inside
    ParaView, while for MEDCoupliung and HDF5 file formats you need additional
    plugin readers.

    There is no plan to make this class complex, it should be as simple as
    possible.

    """
    def __init__(self, file_path: str=""):
        self.mesh_name: str = "mesh"

        self.vertices: np.ndarray = np.array([])
        self.triangles: np.ndarray = np.array([])

        # Field index
        self.index = 0
        self.time = 0.0
        self.number_of_time_steps = 1

        # Fields
        self.arrays: Dict[str, Dict[int, np.ndarray]] = {}

        # Writing data
        self.info_on_components: Dict[str, Dict[int, list]] = {}
        self.times: Dict[str, Dict[int, float]] = {}

        # MODULE pointer
        self.backend = None
        self.backend_name = 'UNKNOWN'
        if file_path:
            self.file_path = file_path
            self.openFile(self.file_path)
        else:
            self.file_path = ""



    def openFile(self, file_path: str):
        """This function does not really open the file because depending on the
        backend, you can either read/write at once, meaning you do not have a
        file handle or you can have it. The backend submodules handle that.
        """
        # Get the extension
        ext = os.path.splitext(file_path)[-1].lower()
        if ext.startswith("."):
            ext = ext[1:]
        self.file_path = file_path

        if ext == "med":
            from . import _medcoupling
            self.backend = _medcoupling
            self.backend_name = "MED"
        elif ext.startswith('vt'):
            from . import _vtk
            self.backend = _vtk
            self.backend_name = "VTK"
        elif "h5" in ext or "hdf" in ext:
            from . import _hdf5
            self.backend = _hdf5
            self.backend_name = "HDF5"
        else:
            self.backend_name = "UNKNOWN"
            self.backend = None

    def isOpen(self) -> bool:
        if self.backend_name == "UNKNOWN":
            return False
        return True

    def getMeshData(self) -> Tuple[np.ndarray]:
        """Get the vertices and triangles of the mesh.
        """
        if self.vertices.size == 0:
            self.readMeshData()

        return self.vertices, self.triangles

    def readMeshData(self) -> None:
        """Returns 2 arrays, vertices and triangles where the shapes of the
        arrays are (N_vertices, 3) with dtype float32 and (N_triangles, 3) with
        dtype uint32 respectively.
        """
        if self.vertices.size == 0:
            v, t = self.backend.getMeshData(self.file_path)
            self.vertices = v
            self.triangles = t

    def setIndex(self, index: int):
        """Sets the index for reading fields.
        """
        self.index = index

    def setTime(self, time: int):
        """Sets the time for reading fields.
        """
        self.time = time

    def setNumberOfTimeSteps(self, number_of_time_steps: int) -> None:
        """Define how many time steps are there in the file. Use this when
        writing data to a file. Some writers (VTK) requires this information
        before writing the data.

        Arguments:
            number_of_time_steps (int): Number of time steps.
        """
        if number_of_time_steps <= 0:
            return
        self.number_of_time_steps = number_of_time_steps

    def getNumberOfTimeSteps(self) -> int:
        """Get number of time steps in a FLT mesh file. It will look for the
        most common FLT field, which is conlen.

        Returns:
            number_of_fields(int): Number of time fields in the file.
        """
        if self.backend is None:
            return 0
        flt_field = "conlen"
        return self.backend.getNumberOfTimeSteps(self.file_path, flt_field)

    def getField(self, field_name, index=None) -> np.ndarray:
        """Returns a field of the the opened file. Index is always provided,
        where it signifies the time step of the field.

        Returns:
            field (np.ndarray): A numpy array of the field.
        """
        if index is None:
            index = self.index
        if self.backend:
            return self.backend.getField(self.file_path, field_name, index)

    def writeMeshTo(self, file_path: str) -> 'Mesh':
        """In this case we wish to copy the mesh of a file to a new location.
        Same file paths are *not* allowed, i.e., over-writing mesh in the same
        file.

        The function returns a new Mesh object that has a copy of only mesh
        data.

        Arguments:
            file_path (str): Path to new file.

        Returns:
            mesh_obj (Mesh): Mesh object with reference to the new file path.

        """
        if file_path == self.file_path:
            return None

        if not self.vertices.size:
            # Try getting the data
            self.readMeshData()
            if not self.vertices.size:
                return None

        # Figure out the correct backend.
        ext = os.path.splitext(file_path)[-1].lower()
        if ext.startswith("."):
            ext = ext[1:]
        backend = None
        if ext == "med":
            from . import _medcoupling
            backend = _medcoupling
        elif ext == "vtu":
            from . import _vtk
            backend = _vtk

        if backend is None:
            return None
        backend.writeMesh(file_path, self.vertices, self.triangles,
                          self.mesh_name)

        # Make a copy
        new_obj = Mesh()
        new_obj.mesh_name = self.mesh_name
        # Copy only mesh
        new_obj.vertices = np.copy(self.vertices)
        new_obj.triangles = np.copy(self.triangles)
        new_obj.openFile(file_path)
        return new_obj

    def addField(self, array_name: str, array: np.ndarray,
                 info_on_components: list=[]):
        """Add field to be written to a file. The input array is actually a
        numpy array, but when we are talking on saving it to files, where the
        array corresponds to an unstructured grid made of triangles, it is
        named a field.

        Arguments:
            array_name (str): Name of the field. It is named array as the data
                is an array, it becomes a field when written in the file.
            array (np.ndarray): Numpy array. Can be multi-dimensional, i.e.,
                vector data.
            index (int): Time step of the array.
            time (float): Self-describing.
            info_on_components (list): When we have vector data it is necessary to also
                provides the info_on_components of each component of the vector data.
        """
        if array_name not in self.arrays:
            self.arrays[array_name] = {}

        if array_name not in self.info_on_components:
            self.info_on_components[array_name] = {}

        if array_name not in self.times:
            self.times[array_name] = {}
        index = self.index
        self.arrays[array_name][index] = array
        self.info_on_components[array_name][index] = info_on_components
        self.times[array_name][index] = self.time

    def writeFields(self):
        """Write the fields to the file.
        """
        if self.backend_name == "MED":
            for array_name in self.arrays:
                for index in self.arrays[array_name]:
                    self.backend.writeField(self.file_path, array_name,
                        self.arrays[array_name][index], index,
                        self.times[array_name][index],
                        self.info_on_components[array_name][index])
        elif self.backend_name == "VTK":
            # Haven't found a way to append array data to an XML Unstructured
            # Grid data, though if I wrote my own XML writer it would be
            # easier. In this case we only the VTU, or the XML Unstructured
            # Grid file is supported as it can actually contains all time steps
            # in a single file. Other formats requires additional XML, or files
            # to be generated.


            # Manually do the thing as it is quite compilicated.
            writer = self.backend.getWriter(self.file_path)
            vtk_orig_obj = self.backend.getVtkObj(self.file_path, self.vertices,
                                                 self.triangles)

            if self.index == 0:
                # We are starting!
                writer.SetFileName(self.file_path)
                writer.SetNumberOfTimeSteps(self.number_of_time_steps)
                writer.Start()

            vtk_obj = self.backend.copyVtkObj(vtk_orig_obj)

            # The problematic code is following, where we have to dump the
            # arrays to the VTK object.
            for array_name in self.arrays:
                self.backend.addArrayToVtkObj(vtk_obj, array_name,
                    self.arrays[array_name][self.index])

            # Now write
            writer.SetInputData(vtk_obj)
            writer.WriteNextTime(self.time)

            if self.index == self.number_of_time_steps - 1:
                writer.Stop()


        return True

    def translateMesh(self, vector: np.ndarray) -> None:
        """Translates the vertices with the provided vector.

        Arguments:
            vector (np.ndarray): Direction vector in which the mesh points are
                translated.
        """
        if self.file_path is None:
            log.error("Could not read mesh data as no MED file is loaded")
            return None

        if self.vertices.size == 0:
            self.readMeshData()
        self.vertices += vector

    def rotateMesh(self, p1: np.ndarray, p2: np.ndarray, angle: float) -> 'Mesh':
        """Rotates the the mesh data around the provided axes and angle. This
        rotation is not "free form", the rotational axis is set to the
        geometry's center of mass.


        It creates a new Mesh, since it creates a copy in any case and it
        is sometimes useful if we retain the original data.

        Arguments:
            p1 (np.ndarray): First point of the rotational vector axis
            p2 (np.ndarray): Second point of the rotational vector axis
            angle (float): Angle of rotationa

        Returns:
            mesh (Mesh): New mesh object with rotated vertices.
        """

        if self.file_path is None:
            log.error("Could not read mesh data as no MED file is loaded")
            return None

        if self.vertices.size == 0:
            self.readMeshData()

        out = Mesh()
        vertices = self.vertices.copy()

        # Now rotate
        new_vertices = rotatePointsAroundAxis(vertices, p1, p2, angle)
        out.vertices = new_vertices
        out.triangles = self.triangles.copy()
        return out

def load_flt_results_from_mesh(mesh_results: L2GResults, mesh_object: Mesh):
    """Loads required fields from the Mesh objects to the mesh results.

    It is necessary to set the index of the fields *before* calling this \
    function.
    """
    if not mesh_object.isOpen():
        log.error("Mesh file is not open to load data from")
        return

    mesh_results.drsep = mesh_object.getField('drsep') * 1e-3 # Convert to m
    mesh_results.conlen = mesh_object.getField('conlen')
    mesh_results.flux = mesh_object.getField('flux')

    # Necessary for the mask.
    mesh_results.geom_hit_ids = mesh_object.getField("geom_hit_ids")

    # Necessary for FL
    # Necessary to flatten!
    mesh_results.BVec = mesh_object.getField("BVec").flatten()
    mesh_results.Bdot = mesh_object.getField("Bdot")
    # Necessary to flatten!
    mesh_results.baryCent = mesh_object.getField("baryCent").flatten()

def dump_flt_results_to_mesh(mesh_results: L2GResults, mesh_object: Mesh):
    """Dumps flt_mesh results to mesh object.

    mesh_results is L2GResults

    """

    for key in mesh_results.arrays_to_dump:
        array = getattr(mesh_results, key)
        if array is None:
            continue

        if key == "angle":
            # Transform the angle
            array = np.rad2deg(array)
            array = np.where(array > 90.0, array - 90.0, 90.0 - array)

        if key in ["drsep", "drsep2", "conlen"]:
            # Scale to mm from m
            array = array * 1000

        infoOnComponent = []
        if key == 'BVec':
            # Write the Magnetic field vector properly. That is a
            # 3 component array in Cartesian!
            infoOnComponent = ['x', 'y', 'z']
        elif key == 'normals':
            infoOnComponent = ['x', 'y', 'z']
        elif key == 'BVecCyln':
            # Write the Magnetic field poloidal and toroidal component
            # properly. A 2 component array in Cylindrical
            infoOnComponent = ['Pol', 'Tor']
        elif key == 'baryCent':
            # Again, write the baryCent properly. In this case the
            # infoOnComponent will basically tell the user that the vector
            # is in Cylindrical system.
            infoOnComponent=['r', 'z', 'phi']

        mesh_object.addField(array=array, array_name=key,
                             info_on_components=infoOnComponent)

    mesh_object.writeFields()

def dump_hlm_results_to_mesh(hlm_results: L2GResultsHLM, mesh_object: Mesh,
                             hlm_type: str):
    """Dumps flt_mesh results to mesh object.

    hlm_results is L2GResultsHLM

    """
    # Classical arrays.
    mesh_object.addField(array=hlm_results.flux_expansion,
                         array_name="Total_flux expansion")
    mesh_object.addField(array=hlm_results.q_inc,
        array_name=r"$q_{\perp}\;\;[\frac{W}{m^2}]$")
    mesh_object.addField(array=hlm_results.q_par,
        array_name=r"$q_{\parallel}\;\;[\frac{W}{m^2}]$")

    if hlm_type == "elm":
        # Writing additional arrays
        # Take arrays from additional_arrays variable.
        mesh_object.addField(
            array=hlm_results.additional_arrays[0],
            array_name=r"$q_{\parallel}\;\;[\frac{W}{m^2}]$ ELM")
        mesh_object.addField(
            array=hlm_results.additional_arrays[1],
            array_name=r"$q_{\parallel}\;\;[\frac{W}{m^2}]$ inter-ELM")
        mesh_object.addField(
            array=hlm_results.additional_arrays[2],
            array_name=r"$T_e$ [eV]")
        mesh_object.addField(
            array=hlm_results.additional_arrays[3],
            array_name=r"$T_i$ [eV]")
        mesh_object.addField(
            array=hlm_results.additional_arrays[4],
            array_name=r"$q_{\perp}\;\;[\frac{W}{m^2}]$ ELM")
        mesh_object.addField(
            array=hlm_results.additional_arrays[5],
            array_name=r"$q_{\perp}\;\;[\frac{W}{m^2}]$ inter-ELM")

    mesh_object.writeFields()

def save_fls_to_vtk(data, file_path):
    """Helper function that saves field lines points or polylines to a VTK
    file. In this case ``vtkGenericDataObjectWriter`` is used.

    This should be used to save :class:`L2GFLs` to a VTK file.
    """
    import vtk
    writer = vtk.vtkGenericDataObjectWriter()
    writer.SetInputData(data)
    writer.SetFileName(file_path)
    writer.Write()

def save_mesh_to_vtk(data, file_path):
    """Helper function that saves field lines points or vtkUnstructuredGrid to
    a VTK file. In this case ``vtkXMLUnstructuredGridWriter`` is used.

    This should be used to save :class:`L2GResults` to a VTK file.
    """
    import vtk
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(data)
    writer.SetFileName(file_path)
    writer.Write()

def save_results_to_vtk(result_obj: Union[L2GFLs, L2GResults],
                       file_path: str) -> None:
    """Helper function to save result objects: :class:`L2GResults`,
    :class:`L2GFLs` to a file.

    In case of :class:`L2GResults` the data is saved to an existing MED file,
    so the MED file has to exist.

    In case of :class:`L2GFLs` a new VTK file is created

    """
    if not isinstance(result_obj, (L2GResults, L2GFLs)):
        log.error("Wrong object provided to results")

    data = result_obj.generate_vtk_object()

    if isinstance(result_obj, L2GResults):
        save_mesh_to_vtk(data, file_path)
    elif isinstance(result_obj, L2GFLs):
        save_fls_to_vtk(data, file_path)