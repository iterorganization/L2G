"""This module contains classes and functions used by the FLT for retrieving
and writing mesh data and/or results.
"""
import numpy as np
import l2g.utils.meshio
from l2g.comp import L2GResults, L2GFLs
import medcoupling as mc
import os

import logging
log = logging.getLogger(__name__)

from typing import Union

class MEDMeshIO():
    """Wrapper class for IO of MED object/files. Also to hold the
    information about the file path and reference to MED mesh object.

    To correctly apply the index and associated time of a field in the mesh,
    users have to edit the attributes :attr:`index` and :attr:`associated_time`

    Warning: Currently we have one mesh per file. In future we might support
    different handling of meshes inside a MED file.
    """

    def __init__(self):
        self.med_mesh = None
        self.med_file_orig_path = "" # Usually we make a copy of input mesh.
        self.med_file_path = "" # Target location of the resulting mesh.

        self.overwrite_mesh = True

        # Used to associate the index and associated time to a slice of data.
        self.index = 0 # Default
        self.associated_time = -1 # Default

        self._vertices = None
        self._triangles = None

    def loadMedFile(self, file: str) -> None:
        if not os.path.exists(file):
            log.error(f"File {file} does not exist!")
            return
        self.med_mesh = mc.ReadMeshFromFile(file)
        self.med_file_orig_path = file

    def getMeshData(self):
        """Return the vertices and triangles of the mesh file.

        Not implemented: Retrieve the data from a submesh of MED file.
        """
        if self.med_mesh is None:
            log.error("Could not read mesh data as no MED file is loaded")
            return None, None

        vertex = self.med_mesh.getCoords().toNumPyArray()
        vertex = np.array(vertex, dtype=np.float32)

        # Triangles

        Ntriangles = self.med_mesh.getNumberOfCellsWithType(mc.NORM_TRI3)

        triangles = np.empty(Ntriangles * 3, np.uint32)

        Ncells = self.med_mesh.getNumberOfCells()

        c = 0
        for i in range(Ncells):
            if self.med_mesh.getTypeOfCell(i) != mc.NORM_TRI3:
                continue

            triangles[c:c+3] = self.med_mesh.getNodeIdsOfCell(i)
            c += 3

        return vertex.reshape(vertex.shape[0] * 3), triangles

    def writeMesh(self, file_location: str) -> None:
        """Write the mesh into a new location. Reason behind is that usually
        we have input mesh and the resulting mesh will be a copy of that mesh.
        """
        if self.med_mesh is None:
            log.error("Could not write mesh as no MED mesh was read!")
            return

        self.med_file_path = file_location
        mc.WriteMesh(fileName=self.med_file_path, mesh=self.med_mesh,
                     writeFromScratch=self.overwrite_mesh)

    def writeArray(self, array: np.ndarray, array_name: str,
                   info_on_component: list = []):
        """Helper function to write a numpy array to the mesh, as several steps
        need to be performed before writing to the med file.
        """
        if self.med_mesh is None:
            log.error(f"Could not write {array_name} as no MED mesh was read!")
            return

        field = l2g.utils.meshio.numpyArrayToField(arr=array,
            field_name=array_name,
            mesh=self.med_mesh, iteration=self.index,
            associated_time=self.associated_time,
            info_on_component=info_on_component)
        l2g.utils.meshio.writeFieldToAlreadyExistingMesh(field=field,
            med_file=self.med_file_path)

    def getArray(self, field_name):
        """Helper function for obtaining numpy array from a MED field.
        """
        return l2g.utils.meshio.fieldToNumpy(self.med_file_path,
                                             field_name, self.index, -1)

def load_flt_mesh_results_from_med(mesh_results: L2GResults, med_object: MEDMeshIO):
    """Loads the FLT data from med object. Usually we only need a few FLT
    relevant quantities.
    """
    if not med_object.med_file_path:
        log.error("No MED file path to load data from")
        return
    mesh_results.drsep = med_object.getArray("drsep")
    mesh_results.conlen = med_object.getArray("conlen")

def dump_flt_mesh_results_to_med(mesh_results: L2GResults, med_object: MEDMeshIO):
    """Dumps flt_mesh results to med object.

    Index and associated time has to be set by user before hand.

    """

    for key in mesh_results.__slots__:
        array = getattr(mesh_results, key)
        if array is None:
            continue

        if key == "angle":
            # Transform the angle
            array = np.rad2deg(array)
            array = np.where(array > 90.0, array - 90.0, 90.0 - array)

        if key in ["drsep", "drsep2", "conlen"]:
            # Scale to mm from m
            array *= 1000

        infoOnComponent = None
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

        med_object.writeArray(array=array, array_name=key,
                              info_on_component=infoOnComponent)

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