import h5py
import numpy as np
from typing import Dict

"""

# Create the HDF5 file
with h5py.File("output.h5", "w") as file:
    # Create the unstructured grid group
    grid_group = file.create_group("unstructured_grid")

    # Create and save the vertices and triangles of the mesh
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], ...])
    triangles = np.array([[0, 1, 2], [1, 2, 3], ...])

    grid_group.create_dataset("vertices", data=vertices)
    grid_group.create_dataset("triangles", data=triangles)

    # Create the time-varying data group
    time_group = file.create_group("time_data")

    # Loop over each time step
    for time_step in range(num_time_steps):
        # Create the datasets for the time-varying fields
        time_step_group = time_group.create_group(f"time_step_{time_step}")

        # Create and save the attribute data for the current time step
        attribute_data = np.random.rand(num_cells)  # Replace with your own attribute data
        time_step_group.create_dataset("attribute", data=attribute_data)

        # ... Create other datasets for time-varying fields ...

        # Store the time value as an attribute
        time_step_group.attrs["time_value"] = time_step * delta_t  # Adjust as needed
"""

class HDF5Backend(object):
    """Backend class for reading/writing unstructured triangular mesh and
    time-varying field data linked on that file.
    """
    def __init__(self, file_path: str=""):
        self.time: int = 0
        self.index: int = 0

        self.backend_obj = None
        self.error = False

        # Fields
        self.arrays: Dict[str, Dict[int, np.ndarray]] = {}

        # Writing data
        self.units: Dict[str, Dict[int, list]] = {}
        self.times: Dict[str, Dict[int, float]] = {}

        if file_path:
            self.file_path: str = file_path
            self.openFile(file_path)
        else:
            self.file_path: str = ""

    def __dell__(self):
        if self.backend_obj:
            self.backend_obj.close()

    def openFile(self, file_path: str):
        """Open a handle to the file.
        """
        self.backend_obj = h5py.File(file_path)
        pass

    def isOpen(self):
        return True

    def checkIfMeshDataIsPresent(self):
        if not self.isOpen():
            return False

        # Check if unstructured_grid is inside
        if "unstructured_grid" in self.backend_obj:
            ugrid = self.backend_obj["unstructured_grid"]

            if "vertices" not in ugrid:
                return False
            if "triangles" not in ugrid:
                return False
        else:
            return False

    def getMeshData(self):
        """Get mesh data of the open file.

        It returns two 1D arrays of vertices and triangles respectively.
        """

        if self.checkIfMeshDataIsPresent():
            return [], []

        ugrid = self.backend_obj["unstructured_grid"]
        # Reshape into 1D array
        vertices = np.asarray(ugrid["vertices"], dtype=np.float32).flatten()

        triangles = np.asarray(ugrid["triangles"], dtype=np.uint32).flatten()
        return vertices, triangles

    def getArray(self, array_name: str, index: int):
        """Function for obtaining a quantity. In the kwargs it really
        """
        if not self.isOpen():
            return []

        # Check for data
        if not "data" in self.backend_obj:
            return []
        data = self.backend_obj["data"]
        time = str(int(index)) # Just to make sure, we convert it to int again

        if time not in data:
            return []

        time_data = data[time]

        if array_name not in time_data:
            return []

        array = np.asarray(time_data[array_name], np.float32)

        return array

    def writeMesh(self, file_path: str, vertices: np.ndarray,
                  triangles: np.ndarray):
        backend_closed = False
        if self.file_path == file_path and self.isOpen():
            # Maybe it's good that we close the already open file
            if self.checkIfMeshDataIsPresent():
                # Data is there and overwrite is not
                return

            # Let's close the backend_obj as it is in R mode.
            self.backend_obj.close()
            backend_closed = True

        with h5py.File(file_path, "w") as f:
            grid_group = f.create_group("unstructured_grid")
            grid_group.create_dataset("vertices", data=vertices)
            grid_group.create_dataset("triangles", data=triangles)

        if backend_closed:
            self.backend_obj = h5py.File(file_path, "r")

    def addArray(self, array_name: str, array: np.ndarray, index:int = 0,
                 time: float = 0.0, units: list = []):
        if array_name not in self.arrays:
            self.arrays[array_name] = {}

        if array_name not in self.units:
            self.units[array_name] = {}

        if array_name not in self.times:
            self.times[array_name] = {}

        self.arrays[array_name][index] = array
        self.units[array_name][index] = units
        self.times[array_name][index] = time

    def writeArray(self, **kwargs):
        """Write a field of one quantity to the file.
        """

        pass
