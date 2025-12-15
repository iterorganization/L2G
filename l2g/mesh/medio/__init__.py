import logging
import os
import h5py
import numpy as np

from .write import (create_group_maps, create_field_structure,
                    create_basic_mesh_structure, add_med_info, write_mesh_data,
                    create_mesh_group_entry, write_mesh_fam_id,
                    write_mesh_group_data, check_if_field_compatible,
                    create_field_entry, write_field_one_step_data)
from .read import (traverse_med_file, get_field, get_mesh_data,
                   get_mesh_fam_id_array, Mesh)

log = logging.getLogger(__name__)

class MEDTR3Reader(object):
    """Python class for reading MED file created by medfile library. More
    specifically it is meant for reading the mesh data, groups and fields of
    triangular surface meshes.

    .. code-block:: python

       import l2g.mesh.medio
       reader = l2g.mesh.medio.MEDTR3Reader('/path/to/hdf5/med/file')
       ...

    """
    def __init__(self, file_path: str):
        """Provide a path to a MED hdf5 file.

        Arguments:
            file_path (str): Path to MED hdf5 file.
        """
        # Create dictionaries, as when we
        self.meshes: dict[str, Mesh] = {}
        self.mesh_fam_ids: dict[str, np.ndarray] = {}
        self.fields: dict[str, list[tuple[int, float]]] = {}
        self.mesh_names: list[str] = []
        self.file_path: str = file_path
        if not os.access(file_path, os.R_OK):
            raise Exception(f"Can't access file {file_path}")

        self.traverse()

    def traverse(self) -> None:
        """Using h5py traverse the MED file to figure out how many meshes and
        groups are there. This function is run on object instance.

        The logic of the function is on-hand experience.

        Right now it obtains information of all meshes and the groups from the
        MED file.

        Meshes:
        A single MED file can have multiple meshes stored with groups. This
        function gathers the names of meshes and the groups tied to it.

        Groups:
        When finding out group names, the name is hidden inside
        FAS/mesh_name/indexed_group_name/GRO/NOM. The NOM is a int array of
        shape ((N,80)) of ascii codes. Integer 0 acts as null terminator for
        the string. If N > 1 it means the following:
         - First name is the name of the group
         - Second and on names are the names of other groups that contain this
           group. At least that' the theory.

        """

        with h5py.File(self.file_path, 'r') as f:
            mesh_names, meshes, fields = traverse_med_file(f)

            self.mesh_names = mesh_names
            self.meshes = meshes
            self.fields = fields

    def getAllMeshes(self) -> list[str]:
        """Return the names of the meshes inside a MED file. Use this to
        determine if there are multiple meshes in a file.

        Returns:
            meshes (list[str])

        """
        out = []

        for mesh_name in self.meshes:
            out.append((mesh_name, list(self.meshes[mesh_name].groups.keys())))

        return out

    def getAllFields(self, mesh_name: str = "") -> list[str]:
        """Returns all fields found inside the MED file.

        Returns:
            fields (list[str])
        """
        if mesh_name == "" and len(self.meshes) > 1:
            log.debug("Warning: Mesh name not provided and MED file contains more than one mesh!")
            # Python >3.6 order of inserted keys in dictionary is preserved
            mesh_name = list(self.meshes.keys())[0]
            log.debug(f"Taking the first occurring mesh: {mesh_name}")
        elif mesh_name != "" and mesh_name not in self.meshes:
            raise Exception(f"Mesh with the name {mesh_name} not in MED file {self.file_path}!")
        else:
            mesh_name = self.mesh_names[0]

        return list(self.fields.keys())

    def getAllFieldIterations(self, field_name: str) -> list[tuple[int, float]]:
        """Returns the [indexes], [times] arrays of the field.

        Returns:
            index (list[int]): List of integers
            times (list[float]): List of floats corresponding to time of the
                same index.
        """
        if field_name not in self.fields:
            raise Exception(f"No field {field_name} in MED file {self.file_path}")
        return self.fields[field_name]

    def getField(self, field_name: str, index: int, mesh_name: str = "") -> np.ndarray:
        """Returns a field associated to a mesh and index from the MED file.

        Arguments:
            field_name (str): Name of a field
            index (int): Index of a field
            mesh_name (str): Name of mesh. Optional, if there is only one mesh
                in a MED file, leave an empty string.

        Returns:
            field (np.ndarray): A (dim, N) numpy array. Dim corresponds to the
                NCO attribute and N corresponds to the NBR attribute
        """
        if field_name not in self.fields:
            raise Exception(f"Field {field_name} not in MED file {self.file_path}")

        if mesh_name == "" and len(self.meshes) > 1:
            log.debug("Warning: Mesh name not provided and MED file contains more than one mesh!")
            # Python >3.6 order of inserted keys in dictionary is preserved
            mesh_name = list(self.meshes.keys())[0]
            log.debug(f"Taking the first occurring mesh: {mesh_name}")
        elif mesh_name != "" and mesh_name not in self.meshes:
            raise Exception(f"Mesh with the name {mesh_name} not in MED file {self.file_path}!")
        else:
            mesh_name = self.mesh_names[0]

        # Field is here, but let's check if there is a mesh with this field
        if not field_name in self.meshes[mesh_name].fields:
            raise Exception(f"Field {field_name} is not tied to mesh {mesh_name}")

        with h5py.File(self.file_path, "r") as f:
            field = get_field(f, field_name, mesh_name, index)
        return field

    def getMeshData(self, mesh_name: str = "") -> tuple[np.ndarray, np.ndarray]:
        """Returns the mesh data, notably the coordinates and the triangles
        indices.

        MEDCoupling by itself takes a lot of parameters, such as

        What is assumed:
            - The MED file has a triangular unstructured mesh inside
            - The Data is accessed by manually hard-coding parts of the path,
              such as injecting "-0000000000000000001-0000000000000000001" into
              the path to obtain mesh data as I haven't seen much usage of
              these two parameters (integers) in practical examples.
            - Data is stored in row Major order. So an array of 3 dimensions
              have the data stored as A = [x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3]
              and to change it into column major, reshape it as:
                                    reshape((3, NBR)).T

        Arguments:
            mesh_name (str): Name of mesh. Optional, if there is only one mesh
                in a MED file, leave an empty string.

        Returns:
            coords (np.ndarray): Numpy array of shape (n coords, 3)
            vertices (np.ndarray): 1D numpy array

        """

        if mesh_name == "" and len(self.meshes) > 1:
            log.debug("Warning: Mesh name not provided and MED file contains more than one mesh!")
            # Python >3.6 order of inserted keys in dictionary is preserved
            mesh_name = list(self.meshes.keys())[0]
            log.debug(f"Taking the first occurring mesh: {mesh_name}")
        elif mesh_name != "" and mesh_name not in self.meshes:
            raise Exception(f"Mesh with the name {mesh_name} not in MED file {self.file_path}!")
        else:
            mesh_name = self.mesh_names[0]

        # Obtain the coordinates
        with h5py.File(self.file_path, 'r') as f:

            vertices, triangles = get_mesh_data(f, mesh_name)

        # Triangles are in Fotrtran notation
        return vertices, triangles - 1

    def getGroupArray(self, group_name: str, mesh_name: str = ""):
        """Obtains the element IDs that belong to a group in a mesh inside the
        MED file. Used mainly for masking quantity data on a mesh.

        .. note::

           It's up to the user to know what group are they accessing, i.e., the
           role of the group. Whether those are edges or faces or etc...

        Arguments:
            group_name (str): Name of group
            mesh_name (str): Name of mesh. Optional, if there is only one mesh
                in a MED file, leave an empty string.

        Returns:
            ids (np.ndarray): 1D numpy array of element IDs.
        """

        if mesh_name == "" and len(self.meshes) > 1:
            log.debug("Warning: Mesh name not provided and MED file contains more than one mesh!")
            # Python >3.6 order of inserted keys in dictionary is preserved
            mesh_name = list(self.meshes.keys())[0]
            log.debug(f"Taking the first occurring mesh: {mesh_name}")
        elif mesh_name != "" and mesh_name not in self.meshes:
            raise Exception(f"Mesh with the name {mesh_name} not in MED file {self.file_path}!")
        else:
            mesh_name = self.mesh_names[0]

        mesh = self.meshes[mesh_name]
        if group_name not in mesh.groups:
            raise Exception(f"Mesh {mesh_name} does not contain group {group_name}")


        group_ids = mesh.groups[group_name]

        if mesh_name not in self.mesh_fam_ids:
            with h5py.File(self.file_path, 'r') as f:
                fam = get_mesh_fam_id_array(f, mesh_name)

            self.mesh_fam_ids[mesh_name] = fam
        fam = self.mesh_fam_ids[mesh_name]
        mask = np.zeros(fam.shape)
        for gid in group_ids:
            mask = np.logical_or(mask, fam == gid)

        return mask

class MEDTR3Writer(object):
    """Python class for writing groups, fields and mesh to a HDF5 .med file,
    which aims to be compatible with the MED file library. It is meant for
    writing 3D surface triangular mesh and its groups and fields.

    .. code-block:: python

       import l2g.mesh.medio
       writer = l2g.mesh.medio.MEDTR3Writer('mesh_name')

       # Add mesh data
       # Add groups
       # Add fields

       writer.write('/path/to/file')
       ...

    """
    def __init__(self, mesh_name: str):
        self.mesh_name: str = mesh_name

        self.vertices: np.ndarray = np.ndarray([])
        self.edges: np.ndarray = np.ndarray([])
        self.triangles: np.ndarray = np.ndarray([])

        ## Groups
        self.groups: dict[str, np.ndarray] = {}
        self.groups_processed: bool = False

        self.group_data: list = []
        self.group_fam_mask: np.ndarray = np.ndarray([])

        ## Fields
        self.fields: dict[str, dict[int, tuple[float, np.ndarray]]] = {}
        self.field_components: dict[str, list[str]] = {}
        self.field_names: set = set()

        self.maj_ver: int = 4
        self.min_ver: int = 1
        self.rel_ver: int = 1

    def addMeshData(self, vertices: np.ndarray, edges: np.ndarray,
                    triangles: np.ndarray):
        self.vertices = vertices
        self.edges = edges
        self.triangles = triangles

    def addGroup(self, group_name: str, group_id: np.ndarray):
        self.groups[group_name] = group_id
        self.groups_processed = False

    def registerField(self, field_name: str, components: list[str] = [""]):
        """Register a field by specifying a name and components. Components
        define if a field is vector field, tensor, scalar, etc... Existing
        registered field gets overwritten.

        You must first register a field before calling write, if you added
        field data related to this field.

        Arguments:
            field_name (str): Name of field
            components (list[str]): Components of the field. If left empty, it
                is a scalar field.
        """

        self.field_components[field_name] = components
        self.field_names.add(field_name)

    def addFieldData(self, field_name: str, index: int, time: float,
                     array:np.ndarray):
        """Adds a field data, associated with its index. Time is of descriptive
        nature and array is the contents of the field at that index.

        Arguments:
            field_name (str): Name of field the data corresponds to
            index (int): Index of the data
            time (float): Descriptive information of this particular step.
                Usually described in temporal space, but it can signify values
                of other space.
        """
        if field_name not in self.fields:
            self.fields[field_name] = {}

        self.fields[field_name][index] = (time, array)

    def setVersion(self, maj: int, min: int, rel: int):
        self.maj_ver = maj
        self.min_ver = min
        self.rel_ver = rel

    def write(self, path: str, overwrite: bool=False):
        """Writes mesh, group and field data to a MED file. If file exists data
        gets appended.

        Arguments:
            path (str): Path to file
            overwrite (bool): Overwite field data.
        """

        if not os.access(path, os.F_OK):
            mode = 'w'
        else:
            mode = 'r+'

        if not self.groups_processed:
            gd, gfm = create_group_maps(self.groups, (self.triangles.shape[0],))
            self.group_data = gd
            self.group_fam_mask = gfm

        mesh_name = self.mesh_name

        with h5py.File(path, mode) as f:
            create_basic_mesh_structure(f)
            add_med_info(f, self.maj_ver, self.min_ver, self.rel_ver)
            write_mesh_data(f, mesh_name, self.vertices, self.edges,
                            self.triangles)

            create_mesh_group_entry(f, mesh_name)
            write_mesh_fam_id(f, mesh_name, self.group_fam_mask)

            for (gname, gid, names) in self.group_data:
                write_mesh_group_data(f, mesh_name, gname, gid, names)

            create_field_structure(f)

            # Iterate over fields
            for field_name in self.field_names:
                if field_name not in self.fields:
                    continue

                components = self.field_components[field_name]
                ncomps = len(components)

                create = True
                if overwrite:
                    pass
                elif check_if_field_compatible(f['CHA'], field_name,
                                             mesh_name, ncomps, components):
                    create = False

                if create:
                    create_field_entry(f['CHA'], field_name, mesh_name,
                                       ncomps, components)

                field: h5py.Group = f[f'CHA/{field_name}']

                for index in self.fields[field_name]:
                    time, array = self.fields[field_name][index]
                    write_field_one_step_data(field, array, time, index)

    def writeMesh(self, path: str):
        """Writes mesh data to a MED file.If file exists data gets appended.

        Arguments:
            path (str): Path to MED file.
        """
        if not os.access(path, os.F_OK):
            mode = 'w'
        else:
            mode = 'r+'

        with h5py.File(path, mode) as f:

            create_basic_mesh_structure(f)
            add_med_info(f, self.maj_ver, self.min_ver, self.rel_ver)
            write_mesh_data(f, self.mesh_name, self.vertices, self.edges,
                            self.triangles)

    def writeGroup(self, path: str):
        """Writes group data to MED file. Note that the function does not check
        if the masks are valid or not. Up to the user.

        Arguments:
            path (str): Path to MED file.
        """
        if not os.access(path, os.F_OK):
            mode = 'w'
        else:
            mode = 'r+'

        if not self.groups_processed:
            gd, gfm = create_group_maps(self.groups, (self.triangles.shape[0],))
            self.group_data = gd
            self.group_fam_mask = gfm

        with h5py.File(path, mode) as f:

            create_mesh_group_entry(f, self.mesh_name)
            write_mesh_fam_id(f, self.mesh_name, self.group_fam_mask)

            for (gname, gid, names) in self.group_data:
                write_mesh_group_data(f, self.mesh_name, gname, gid, names)

    def writeField(self, path: str, overwrite: bool=False):
        """Writes field data to a MED file. Destructive behaviour

        Arguments:
            path (str): Path to MED file.
        """
        if not os.access(path, os.F_OK):
            mode = 'w'
        else:
            mode = 'r+'

        mesh_name = self.mesh_name
        with h5py.File(path, mode) as f:

            create_field_structure(f)

            # Iterate over fields
            for field_name in self.field_names:
                if field_name not in self.fields:
                    continue

                components = self.field_components[field_name]
                ncomps = len(components)

                create = True
                if overwrite:
                    pass
                elif check_if_field_compatible(f['CHA'], field_name,
                                             mesh_name, ncomps, components):
                    create = False

                if create:
                    create_field_entry(f['CHA'], field_name, mesh_name,
                                       ncomps, components)

                field: h5py.Group = f[f'CHA/{field_name}']

                for index in self.fields[field_name]:
                    time, array = self.fields[field_name][index]
                    write_field_one_step_data(field, array, time, index)

"""The following functions are sort of a module interface for MED file
operations, even though they do not make the most of sense, as they use the
above defined classes.

One important notice is that in general MED file format supports multiple
meshes inside a single file. However, certain functions below operates only on
the first mesh that was read from the file.
"""

__medtr3reader__: dict[str, MEDTR3Reader] = {}

def getMeshData(file_path: str, mesh_name:str="") -> tuple[np.ndarray, np.ndarray]:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)

    mesh = __medtr3reader__[file_path]
    return mesh.getMeshData()

def getNumberOfTimeSteps(file_path: str, field: str) -> int:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]

    return len(mesh.getAllFieldIterations(field))

def checkIfFieldExists(file_path: str, field: str) -> bool:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    return field in mesh.fields

def getAllFieldNames(file_path: str) -> list[str]:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    return list(mesh.fields.keys())

def getAllFieldIterations(file_path: str, field: str) -> list[tuple[int, float]]:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    return mesh.getAllFieldIterations(field)

def getField(file_path: str, field: str, index: int) -> np.ndarray:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    return mesh.getField(field, index)

def doesItContainGroup(file_path: str, group_name: str) -> bool:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]

    # Only first mesh is checked
    if len(mesh.mesh_names) > 1:
        return group_name in mesh.meshes[mesh.mesh_names[0]].groups
    return False

def getGroupArr(file_path: str, group_name: str) -> np.ndarray:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    return mesh.getGroupArray(group_name)

def getAllGroups(file_path: str) -> list[str]:
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    if len(mesh.mesh_names) == 0:
        return []

    return list(mesh.meshes[mesh.mesh_names[0]].groups.keys())

def writeMesh(file_path: str, vertices: np.ndarray, triangles: np.ndarray,
              mesh_name: str) -> None:

    obj = MEDTR3Writer(mesh_name)
    obj.addMeshData(vertices, np.array([]), triangles)
    obj.write(file_path)

def writeField(file_path: str, field_name: str, array:np.ndarray,
               index: int, time: int, components: list[str]):
    global __medtr3reader__

    if file_path not in __medtr3reader__:
        __medtr3reader__[file_path] = MEDTR3Reader(file_path)
    mesh = __medtr3reader__[file_path]
    mesh_name = mesh.mesh_names[0]

    obj = MEDTR3Writer(mesh_name)
    obj.registerField(field_name, components)
    obj.addFieldData(field_name, index, time, array)
    obj.writeField(file_path)