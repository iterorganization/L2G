import os
import h5py
import numpy as np

class Mesh(object):
    """Since a MED file can contain multiple meshes, a class or struct is used
    to store information and mappings of the groups and other data.

    """
    def __init__(self):
        self.name: str = ""
        # Groups have an ID mask for the elements. See in the traverse function
        # how it is done. If a group does not have an ID (it's has a
        # FAMILLE_ZERO defined then it contains all the elements.)
        self.groups: dict[str, set[int]] = {}
        self.fields: list = []

def traverse_med_file(f: h5py.File) -> tuple[list[str], dict[str, Mesh], dict[str, list[int]]]:
    """Reads a MED HDF5 file and gathers information about meshes, mesh groups
    and fields.

    Arguments:
        f (h5py.File): HDF5 file handle

    Returns:
        mesh_names (list[str]): List of mesh names inside a MED file.
        meshes (dict[str, Mesh]): Dictionary mapping mesh name to mesh structs
            containing information on groups, etc...
        fields (dict[str, list[int]]): Dictionary holding name of fields and
            the iterations (basically what time steps are in MED file.)
    """
    meshes: dict[str, Mesh] = {}
    mesh_names: list[str] = []
    fields: dict[str, list[tuple[int, float]]] = {}

    if (FAS := f.get('FAS')):

        for key in FAS.keys():
            mesh_name = key

            mesh = Mesh()
            mesh.name = mesh_name
            meshes[mesh_name] = mesh

            MESH = FAS[key]
            meshes[mesh_name] = mesh
            mesh_names.append(mesh_name)

            # Processing GROUPS
            if not (ELEME := MESH.get('ELEME')):
                continue

            mesh_groups = list(ELEME.keys())

            # Name in the group is not relevant, probably just a way to create
            # unique name? Actual names are stored in NOM dataset
            for FAM_GROUP in mesh_groups:
                NOM = ELEME[f'{FAM_GROUP}/GRO/NOM']

                # NOM contains the actual labels of the Group. It can
                # contain multiple names, I assume that it means that
                # first name is the name of the group and other names
                # are the names of the group that contains this one

                # Convert the numpy array from shape ((1, 80)) to a
                # 1D array of size 80.

                # NOM_AR = NOM[:]
                ar = list(NOM[:])
                # Convert the array to string. Nullterm is applied here
                #                                                    v
                groups = ["".join(chr(_) for _ in el[:np.where(el==0)[0][0]]) for el in ar]
                # Convert the array from int to ASCII. String is null
                # terminated, so stop at first 0 integer value in ar.
                # groups = "".join([chr(_) for _ in ar[:ar.index(0)]])

                # Group ID is stored in FAM_GROUP.attrs['NUM']
                NUM: int = ELEME[FAM_GROUP].attrs['NUM']
                for group in groups:
                    if group not in mesh.groups:
                        mesh.groups[group] = set()
                    mesh.groups[group].add(NUM)

    # Check for fields
    if (CHA := f.get('CHA')):
        for field_name in CHA.keys():

            iterations: list[tuple[int, float]] = []
            # Get the indexes and times

            # Get the tied mesh name
            mesh_name: str = CHA[field_name].attrs["MAI"].decode()
            meshes[mesh_name].fields.append(str(field_name))

            for entry in CHA[field_name]:
                index = int(entry[:20])
                time = float(CHA[field_name][entry].attrs["PDT"])

                iterations.append((index, time))
            fields[field_name] = iterations

    return mesh_names, meshes, fields

def get_field(f: h5py.File, field_name: str, mesh_name: str, index: int) -> np.ndarray:
    """Get field from a MED file.

    Arguments:
        f(h5py.File): HDF5 handle
        field_name(str): Name of field to obtain
        mesh_name(str): Name of mesh to which it belongs
        index(int): Index of the field.
    """

    # First check if the field exists that is tied to this mesh name
    if f[f'CHA/{field_name}'].attrs['MAI'].decode() != mesh_name:
        raise Exception(f"Field {field_name} tied to mesh {mesh_name} does not exist in {f.filename}")

    # Get number of components
    num_of_components = f[f'CHA/{field_name}'].attrs['NCO']

    # Now check if it exists with the index
    if f"CHA/{field_name}/{index:020d}{-1:020d}" not in f:
        raise Exception(f"Field {field_name} does not have an entry with index {index}")
    path = f"CHA/{field_name}/{index:020d}{-1:020d}/MAI.TR3/MED_NO_PROFILE_INTERNAL/CO"

    number_of_cells = f[f"CHA/{field_name}/{index:020d}{-1:020d}/MAI.TR3/MED_NO_PROFILE_INTERNAL"].attrs['NBR']
    data = f[path][:]

    if num_of_components == 1:
        return data
    else:
        return data.reshape((num_of_components, -1)).T

def get_mesh_data(f: h5py.File, mesh_name: str) -> tuple[np.ndarray, np.ndarray]:
    """Get vertices and triangles of a mesh from MED file.

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
        f(h5py.File): HDF5 file handle.
        mesh_name(str): Name of mesh/

    Returns:
        vertices (np.ndarray): Vertices of mesh of shape ((N, 3))
        triangles (np.ndarray): Triangles of mesh of shape ((N, 3))
    """
    # Obtain the coordinates
    NOE = f[f'ENS_MAA/{mesh_name}/-0000000000000000001-0000000000000000001/NOE/COO']
    # Get the data
    vertices = NOE[:]
    # Reshape
    vertices = vertices.reshape((3, NOE.attrs['NBR'])).T

    # We focus only on TR3 as in triangles and coordinates
    TR3 = f[f'ENS_MAA/{mesh_name}/-0000000000000000001-0000000000000000001/MAI/TR3/NOD']
    triangles = TR3[:]
    triangles = triangles.reshape((3, TR3.attrs['NBR'])).T
    return vertices, triangles

def get_mesh_fam_id_array(f: h5py.File, mesh_name: str) -> np.ndarray:
    """Obtain the family id array (of triangle cells) for mesh.

    Ideally, the name (which is provided), index, order and which cell types
    should also be provided in order to finish this function.

    Arguments:
        f(h5py.File): HDF5 file handle.
        mesh_name(str): Name of mesh/

    Returns:
        fam (np.ndarray): A 1D numpy array containing integers, used to obtain
            group indices.
    """
    FAM = f[f"ENS_MAA/{mesh_name}/-0000000000000000001-0000000000000000001/MAI/TR3/FAM"]
    fam: np.ndarray = FAM[:]
    return fam
