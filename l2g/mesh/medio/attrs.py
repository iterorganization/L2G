"""Module containing HDF5 structures used by MED fichier for writing data to a
MED HDF5 file. Problem is, if you do not follow the same "specialization" then
medfile will return error...
"""

import h5py
import numpy as np


def med_str_dtype() -> object:
    """Get a MED compatible HDF5 String Attribute type. Basically an
    ASCII null-terminated string... MED does not like anything else.

    Argument:
        size (int): Size of string
    """
    t = h5py.h5t.C_S1.copy()
    t.set_strpad(h5py.h5t.STR_NULLTERM)
    t.set_cset(h5py.h5t.CSET_ASCII)
    return t

def create_med_string_attribute(node: h5py.Group | h5py.Dataset, name: str,
                                string: str) -> None:
    """Creates a MED compatible string attribute to a hdf5 node. If same name
    attribute exists, it gets deleted first.

    Arguments:
        node (h5py.Group): HDF5 entry
        name (str): Name of dataset
        string (str): Contents of string
    """
    if name in node.attrs:
        del node.attrs[name]

    data = string.encode('ascii') + b'\x00'
    dtype = med_str_dtype()
    dtype.set_size(len(data))

    # Create scalar space
    sid = h5py.h5s.create(h5py.h5s.SCALAR)
    # Create string attribute
    aid = h5py.h5a.create(node.id, name.encode('ascii'), dtype, sid)
    # Write the data
    aid.write(np.array(data, dtype='S'))
    aid.close()
    sid.close()

def create_med_int32_attribute(node: h5py.Group | h5py.Dataset,
                               name: str, value: int) -> None:
    """Creates a MED compatible int32 attribute to a hdf5 node. If same name
    attribute exists, it gets deleted first.

    Arguments:
        node (h5py.Group): HDF5 entry
        name (str): Name of dataset
        string (str): Contents of string
    """
    if name in node.attrs:
        del node.attrs[name]

    dtype = h5py.h5t.STD_B32LE
    sid = h5py.h5s.create(h5py.h5s.SCALAR)
    aid = h5py.h5a.create(node.id, name.encode('ascii'), dtype, sid)
    aid.write(np.array(value, np.uint32))
    aid.close()
    sid.close()

def create_med_int8le_dataset(node: h5py.Group, name: str, array: np.ndarray) -> h5py.Dataset:
    """Creates a MED compatible int8 (int 8 LE) HDF5 Dataset. If same name
    dataset exists, it gets deleted first.

    This dataset is usually used to store an array of strings converted to
    ASCII code. In other words a list of names of structures or groups.

    Arguments:
        node (h5py.Group): HDF5 entry
        name (str): Name of dataset
        array (np.ndarray): Numpy 1D array
    """
    if name  in node:
        del node[name]
    array_type = h5py.h5t.array_create(h5py.h5t.STD_I8LE, (80,))
    sid = h5py.h5s.create_simple((array.shape[0],))

    did = h5py.h5d.create(node.id, name.encode('ascii'), array_type, sid)
    # Writing the array in the "normal" way doesn't work.
    # did.write(h5py.h5s.ALL, h5py.h5s.ALL, array)
    dataset = h5py.Dataset(did)
    dataset[:] = array
    did.close()
    sid.close()
    return dataset

def create_med_int64le_dataset(node: h5py.Group, name: str, array: np.ndarray) -> h5py.Dataset:
    """Creates a MED compatible int64 (int 64 LE) HDF5 Dataset. If same name
    dataset exists, it gets deleted first.

    Arguments:
        node (h5py.Group): HDF5 entry
        name (str): Name of dataset
        array (np.ndarray): Numpy 1D array
    """
    if name  in node:
        del node[name]

    sid = h5py.h5s.create_simple((array.size,))
    dtype = h5py.h5t.STD_I64LE

    dit = h5py.h5d.create(node.id, name.encode('ascii'), dtype, sid)
    dit.write(h5py.h5s.ALL, h5py.h5s.ALL ,array)
    dit.close()
    sid.close()
    return node[name]

def create_med_float64le_dataset(node: h5py.Group, name: str,
                                 array: np.ndarray) -> h5py.Dataset:
    """Creates a med-file compatible float (float 64 LE) HDF5 Dataset. If
    dataset exists, it gets deleted first.

    Arguments:
        node (h5py.Group): HDF5 entry
        name (str): Name of dataset
        array (np.ndarray): Numpy 1D array
    """

    if name in node:
        del node[name]

    sid = h5py.h5s.create_simple((array.size,))
    dtype = h5py.h5t.IEEE_F64LE
    did = h5py.h5d.create(node.id, name.encode('ascii'), dtype, sid)
    did.write(h5py.h5s.ALL, h5py.h5s.ALL ,array)
    did.close()
    sid.close()
    return node[name]