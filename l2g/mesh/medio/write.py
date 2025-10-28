import numpy as np
import h5py
from .attrs import (create_med_int32_attribute,
                    create_med_string_attribute,
                    create_med_float64le_dataset,
                    create_med_int8le_dataset,
                    create_med_int64le_dataset)

def create_group_maps(groups: dict[str, np.ndarray],
                      shape) -> tuple[list, np.ndarray]:
    """From a given dictionary dict[str, np.ndarray] of boolean masks
    construct an indexed array that maps overlapping of groups if it
    occurs.

    Arguments:
        groups (dict[str, np.ndarray]): A dictionary of either masks or ids
            of a group
        shape (np.ndarray): Shape or size of elements.

    Returns:
        groups (list[tuple[str, int, list[str]]]): List of a group to save into
            MED file. Each entry is HDF5 Group that contains it's unique index
            value and a list of groups that belongs to it.
        mask (np.ndarray): A 1D array of indexes of groups.

    """
    out = []
    MASK = np.zeros(shape)
    SMASK = np.zeros(shape) # Bitmask
    CMASK = np.zeros(shape) # Sum mask
    group_names = list(groups.keys())

    # Construct the mask by bit switching
    bit2name = {}
    for i, name in enumerate(group_names):
        group_mask = groups[name]
        if group_mask.shape == shape:
            SMASK += group_mask << i
            CMASK += group_mask
        else:
            SMASK[group_mask] += 1 << i
            CMASK[group_mask] += 1
        bit2name[i] = name

    famid = -6 # dunno why this

    # First add the groups
    for name in groups:
        _mask = np.zeros(shape)
        # Both a boolean mask or a collection of ids will work.
        if groups[name].shape == shape and groups[name].max() == 1:
            _mask[groups[name].astype(bool)] = 1
        else:
            _mask[groups[name]] = 1
        _mask[CMASK > 1] = 0
        # We have some independent cells
        if np.sum(_mask):
            out.append((f"FAM_{famid}_{name}", famid, [name]))
            MASK[_mask.astype(bool)] = famid
            famid -= 1

    # Now we process the global mask.
    for i in np.unique(SMASK):
        i = int(i)
        if i == 1:
            continue

        idx = np.where(SMASK == i)[0]
        if idx.size:
            # Get the names
            names = [name for bit, name in enumerate(group_names) if (i >> bit) & 1]

            # Not sure we HAVE to create the string as this...
            out.append((f"FAM_{famid}_{'_'.join(names)}", famid, list(names)))
            MASK[idx] = famid
            famid -= 1

    return out, MASK

def create_basic_mesh_structure(f: h5py.File) -> None:
    """Create HDF5 groups at root level that are necessary to store mesh
    information.

    Arguments:
        f(h5py.File): HDF5 file handle.
    """
    keys = f.keys()

    root_elements = ['ENS_MAA', 'FAS', 'INFOS_GENERALES']

    for el in root_elements:
        if el not in keys:
            f.create_group(el)

def add_med_info(f: h5py.File, maj: int=4, min: int=1, rel: int=1):
    """Manipulate the values in order to test if this is the main culprit of
    "incompatibility" between different MED library versions.

    Arguments:
        f (h5py.File): File handle
        maj (int): Major version
        min (int): Minimum version
        rel (int): Rel version
    """
    INFOS_GENERALES = f['INFOS_GENERALES']
    INFOS_GENERALES.attrs.create('MAJ', maj)
    INFOS_GENERALES.attrs.create('MIN', min)
    INFOS_GENERALES.attrs.create('REL', rel)

def create_field_structure(f: h5py.File) -> None:
    """Create Group CHA at root. Track order is required!
    """
    if 'CHA' not in f:
        f.create_group('CHA', track_order=True)

def check_if_field_compatible(f: h5py.Group, field_name: str, mesh_name: str,
                              num_of_components: int = 1,
                              components: list[str] = ['']) -> bool:
    """Checks if an existing HDF5 Group Field data is compatible with what we
    wish to save. This is done by checking the stored reference mesh name,
    number of components the field have and the actual components label. If it
    differs it returns None.

    This function should be used to determine if we wish to delete existing
    Field group, simply because it has incompatible meta data.
    """
    compatible = True
    if field_name in f:
        _attrs = f[field_name].attrs
        if mesh_name not in _attrs or mesh_name != _attrs['MAI']:
            compatible = False

        if 'NCO' not in _attrs or num_of_components != _attrs['NCO']:
            compatible = False

        cstring = "".join(f"{_:<16}" for _ in components)  # pyright: ignore[reportUnusedVariable]
        if 'NOM' not in _attrs or cstring != _attrs['NOM']:
            compatible = False
    else:
        compatible = False
    return compatible

def create_field_entry(f: h5py.Group, field_name: str, mesh_name: str,
                       num_of_components: int = 1,
                       components: list[str] = ['']) -> None:
    """Creates a field entry under /CHA. If the group exists, it gets deleted.
    You only need to create this once. Afterwards you use
    write_field_one_step_data to write actual time step data.

    Arguments:
        f (h5py.Group): HDF5 Node /CHA
        field_name (str): Name of field
        mesh_name (str): Name of mesh the field corresponds to.
        num_of_components (int): Number of components of the field. Higher then
            one means it is a tensor.
        components (list[str]): List of name of components, e.g., ['X', 'Y',
            'Z']

    """
    # Check if it exists. And if it exists check if it's for the same mesh_name

    if field_name in f:
        del f[field_name]

    FIELD = f.create_group(field_name, track_order=True)
    # The LAA and LCA attributes should reflect the number of steps. At least
    # from inspecting MED files, however they are not needed when loading the
    # data in SMESH or opening the MED file in ParaView.
    FIELD.attrs.create('LAA', 1) # Updated after counting the number of actual steps
    FIELD.attrs.create('LCA', 1) # Updated after counting the number of actual steps
    # Some default values?
    create_med_int32_attribute(FIELD, 'LEN', 1)
    create_med_int32_attribute(FIELD, 'LGC', 16)

    FIELD.attrs.create('TYP', 6)

    # Actual information that matters
    FIELD.attrs.create('NCO', num_of_components)
    create_med_string_attribute(FIELD, 'NOM', "".join(f"{_:<16}" for _ in components))
    create_med_string_attribute(FIELD, 'UNI', ' '*16*num_of_components)
    create_med_string_attribute(FIELD, 'UNT', 's')
    create_med_string_attribute(FIELD, 'MAI', mesh_name)

def write_field_one_step_data(f: h5py.Group, values: np.ndarray, dt: float,
    it: int) -> None:
    """Stores data of one time step of a field. If a HDF5 Group entry already
    exists, it gets deleted first.

    Arguments:
        f (h5py.Group): HDF5 Field node
        values (np.ndarray): Values to save
        dt (float): Time
        it (int): Index of the field.
    """

    # See if there is already an entry of the field.
    entry = f"{it:020d}{-1:020d}"

    if entry in f:
        # Remove it
        del f[entry]

    step = f.create_group(entry)

    create_med_int32_attribute(step, 'LEN', 1)
    create_med_int32_attribute(step, 'LGC', 16)
    step.attrs.create('NDT', it)
    step.attrs.create('NOR', -1)
    step.attrs.create('PDT', dt)
    step.attrs.create('RDT', -1)
    step.attrs.create('ROR', -1)

    MAITR3 = step.create_group('MAI.TR3')
    create_med_string_attribute(MAITR3, 'GAU', '')
    create_med_string_attribute(MAITR3, 'PFL', 'MED_NO_PROFILE_INTERNAL')
    MED_NO_PROFILE_INTERNAL = MAITR3.create_group("MED_NO_PROFILE_INTERNAL")
    create_med_string_attribute(MED_NO_PROFILE_INTERNAL, 'GAU', '')
    MED_NO_PROFILE_INTERNAL.attrs.create('NBR', values.shape[0]) # Number of cells
    MED_NO_PROFILE_INTERNAL.attrs.create('NGA', 1)

    create_med_float64le_dataset(MED_NO_PROFILE_INTERNAL, 'CO', values.T.flatten())

def write_mesh_data(f: h5py.File, mesh_name: str, vertices: np.ndarray,
                    edges: np.ndarray, triangles: np.ndarray) -> None:
    """Writes an unstructured grid mesh into a MED file. If a HDF5 Group with
    the same name already exists, it gets deleted.

    Arguments:
        f (h5py.File): A hdf5 file handle
        mesh_name (str): Name of mesh
        vertices (np.ndarray): 2D numpy array of shape (N, 3).
        edges (np.ndarray): 2D numpy array of shape (N, 2).
        triangles (np.ndarray): 2D numpy array of shape (N, 3).
    """

    # Creating mesh
    root: h5py.Group | h5py.Dataset = f['ENS_MAA']
    if mesh_name in root.keys():
        print(f'Deleting existing mesh {mesh_name}')
        del root[mesh_name]

    mesh = root.create_group(mesh_name)

    # Create attributes
    create_med_string_attribute(mesh, 'DES', '')
    mesh.attrs.create('DIM', data= 2)
    mesh.attrs.create('ESP', data= 3)
    create_med_string_attribute(mesh, 'NOM', '')
    mesh.attrs.create('NXI', data= -1)
    mesh.attrs.create('NXT', data= -1)
    mesh.attrs.create('REP', data= 0)
    mesh.attrs.create('SRT', data= 0)
    mesh.attrs.create('TYP', data= 0)
    create_med_string_attribute(mesh, 'UNI', '')
    create_med_string_attribute(mesh, 'UNT', '')

    # Default mesh
    entry = mesh.create_group('-0000000000000000001-0000000000000000001')
    entry.attrs.create('CGT', data= 1)
    entry.attrs.create('NDT', data= -1)
    entry.attrs.create('NOR', data= -1)
    entry.attrs.create('NXI', data= -1)
    entry.attrs.create('NXT', data= -1)
    entry.attrs.create('PDT', data= -1.0)
    entry.attrs.create('PVI', data= -1)
    entry.attrs.create('PVT', data= -1)

    ### Nodes
    NOE = entry.create_group('NOE')
    NOE.attrs.create('CGS', data=1)
    NOE.attrs.create('CGT', data=1)
    create_med_string_attribute(NOE, 'PFL', 'MED_NO_PROFILE_INTERNAL')

    # Coordinates
    COO = create_med_float64le_dataset(NOE, 'COO', vertices.T.flatten())
    COO.attrs.create('CGT', data=1)
    COO.attrs.create('NBR', data=vertices.shape[0])

    # Family
    FAM = create_med_int64le_dataset(NOE, 'FAM', np.zeros(vertices.shape[0]))
    FAM.attrs.create('CGT', data=1)
    FAM.attrs.create('NBR', data=vertices.shape[0])

    # Num. Cell ID starts with 1...
    NUM = create_med_int64le_dataset(NOE, 'NUM', np.arange(vertices.shape[0]) + 1)
    NUM.attrs.create('CGT', data=1)
    NUM.attrs.create('NBR', data=vertices.shape[0])

    ## Mesh
    MAI = entry.create_group('MAI')
    MAI.attrs.create('CGT', data=1)
    if edges.size:
        SE2 = MAI.create_group('SE2')
        SE2.attrs.create('CGS', data=1)
        SE2.attrs.create('CGT', data=1)
        # Edge type?
        SE2.attrs.create('GEO', data=102)
        create_med_string_attribute(SE2, 'PFL', 'MED_NO_PROFILE_INTERNAL')

        # Create FAM
        FAM = create_med_int64le_dataset(SE2, 'FAM', np.zeros(edges.shape[0]))
        FAM.attrs.create('CGT', data=1)
        FAM.attrs.create('NBR', data=edges.shape[0])

        # Create NOD
        NOD = create_med_int64le_dataset(SE2, 'NOD', edges.T.flatten() + 1)
        NOD.attrs.create('CGT', data=1)
        NOD.attrs.create('NBR', data=edges.shape[0])

        # Create NUM
        NUM = create_med_int64le_dataset(SE2, 'NUM', np.arange(edges.shape[0]) + 1)
        NUM.attrs.create('CGT', data=1)
        # Indexing start at 1
        NUM.attrs.create('NBR', data=edges.shape[0])


    if triangles.size:
        TR3 = MAI.create_group('TR3')
        TR3.attrs.create('CGS', data=1)
        TR3.attrs.create('CGT', data=1)
        # Triangle type?
        TR3.attrs.create('GEO', data=203)
        create_med_string_attribute(TR3, 'PFL', 'MED_NO_PROFILE_INTERNAL')

        # Create FAM
        FAM = create_med_int64le_dataset(TR3, 'FAM', np.zeros(triangles.shape[0]))
        FAM.attrs.create('CGT', data=1)
        FAM.attrs.create('NBR', data=triangles.shape[0])

        # Create NOD
        NOD = create_med_int64le_dataset(TR3, 'NOD', triangles.T.flatten() + 1)
        NOD.attrs.create('CGT', data=1)
        NOD.attrs.create('NBR', data=triangles.shape[0])

        # Create NUM
        NUM = create_med_int64le_dataset(TR3, 'NUM', np.arange(triangles.shape[0]) + 1)
        NUM.attrs.create('CGT', data=1)
        # Indexing start at 1
        NUM.attrs.create('NBR', data=triangles.shape[0])

def create_mesh_group_entry(f: h5py.File, mesh_name: str) -> None:
    """Prepares an entry on groups for mesh. Always re-create this data when
    updating/creating groups.

    Arguments:
        f (h5py.File): h5py.File handle
        mesh_name (str): Name of the mesh for the groups.
    """
    FAS: h5py.Group = f['FAS']

    if mesh_name in FAS:
        del FAS[mesh_name]
    mesh = FAS.create_group(mesh_name)

    # Important to have track_order set to true.
    mesh.create_group('ELEME', track_order=True)
    FAMILLE_ZERO = mesh.create_group('FAMILLE_ZERO')
    FAMILLE_ZERO.attrs.create('NUM', 0)

def write_mesh_group_data(f: h5py.File, mesh_name:str, group_name: str,
                          group_index: int, names: list[str]) -> None:
    """Write an entry of group or groups of a mesh. Basically it writes an
    index, which marks elements that belongs to the groups stored in the
    names list. The group_name is nothing, just an entry, at least I do not
    believe it's important the content of group_name string...

    Arguments:
        f (h5py.File): A h5py.File handle for mesh.
        group_name (str): Name of mesh.
        group_index (int): Identifier.
        names (list[str]): Name of lists.
    """

    path = f'FAS/{mesh_name}/ELEME'
    if path not in f:
        # No mesh data exists... Should write that first
        return
    ELEME: h5py.Group = f[path]

    if group_name in ELEME:
        family_group: h5py.Group = ELEME[group_name]
    else:
        family_group = ELEME.create_group(group_name)

    # Create NUM attribute containing the index
    family_group.attrs.create("NUM", group_index)

    # Create an entry for the various groups that share this index
    if 'GRO' in family_group:
        GRO = family_group['GRO']
    else:
        GRO = family_group.create_group('GRO')
    GRO.attrs.create("NBR", len(names))

    # Now create a dataset NOM that is of shape (len(names), 80).

    NOM_array = np.zeros((len(names), 80))
    for i, name in enumerate(names):
        for j in range(min(len(name), 80)):
            NOM_array[i, j] = ord(name[j])

    create_med_int8le_dataset(GRO, "NOM", NOM_array)

def write_mesh_fam_id(file: h5py.File, mesh_name: str,
                      mesh_fam_id: np.ndarray) -> None:
    """Writes the id of groups arrays of triangle cells. See the hardcoded
    MAI/TR3. To 'support' other types, you simply have to create the correct
    group
    """

    path: str = f"ENS_MAA/{mesh_name}/{-1:020d}{-1:020d}/MAI/TR3/FAM"

    if not path in file:
        # Mesh has to be written first in order to make sense to write the
        # elements id data. Maybe raise an exception.
        return

    FAM: h5py.Dataset = file[path]
    FAM[:] = mesh_fam_id
