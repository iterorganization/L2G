import tempfile
import time
import unittest
import os

from l2g.mesh.medio import MEDTR3Reader, MEDTR3Writer
from l2g.mesh import Mesh
import numpy as np

"""Test the l2g.mesh module, by testing the i/o of l2g.mesh.medio and at the
same time testing the l2g.mesh to see if the backend acts normal.
"""

# Mesh data
vertices = np.array([
    [0, 0, 0], # 0
    [1, 0, 0], # 1
    [0, 1, 0], # 2
    [1, 1, 0], # 3
    [2, 0, 0], # 4
    [2, 1, 0], # 5
    [3, 1, 0], # 6
    [3, 0, 0], # 7
    [4, 0, 0], # 8
    [4, 1, 0]  # 9
])
edges = np.array([[]])
triangles = np.array([
    [0, 1, 2], # 0
    [1, 3, 2], # 1
    [1, 4, 3], # 2
    [4, 5, 3], # 3
    [4, 7, 5], # 4
    [7, 6, 5], # 5
    [7, 8, 6], # 6
    [8, 9, 6]  # 7
])
groups = {}

# Mask group
groups["g1"] = np.zeros(triangles.shape[0], dtype=int)
groups["g2"] = np.zeros(triangles.shape[0], dtype=int)
# Idx group
groups["g3"] = np.array([4,5,6,7])

groups["g1"][:4] = 1
groups["g2"][3:7] = 1

x = np.linspace(0, np.pi, triangles.shape[0])
y = np.sin(x)
y2 = np.zeros((triangles.shape[0], 2))
y2[:, 1] == 1

class TestMesh(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        """Create a temporary directory in which to save temporary files.
        """
        cls.tmp_dir = tempfile.TemporaryDirectory()
        cls.test_file = os.path.join(cls.tmp_dir.name, 'test.med')

        cls.triangles = triangles
        cls.vertices = vertices
        cls.groups = groups

        cls.y = y
        cls.y2 = y2

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp_dir.cleanup()

    def test_1_write_mesh(self):
        """Create a test.med file
        """
        writer = MEDTR3Writer('test')
        # Write mesh data and groups
        writer.addMeshData(self.vertices, np.array([]),  self.triangles)
        for group in self.groups:
            writer.addGroup(group, self.groups[group])

        # Now add fields
        writer.registerField('Test')
        writer.addFieldData('Test', 0, 1.5, self.y)
        writer.addFieldData('Test', 1, 3.0, self.y*2)
        writer.addFieldData('Test', 2, 4.22, self.y*3)

        writer.registerField('Test 2', ["1", "2"])
        writer.addFieldData('Test 2', 0, -1, self.y2)
        writer.addFieldData('Test 2', 1, 0, self.y2*2)

        writer.write(self.test_file)
        self.assertTrue(os.access(self.test_file, os.R_OK))

    def test_2_read_mesh_data(self):
        """Read the mesh data from med file created in test 1.
        """
        reader = MEDTR3Reader(self.test_file)

        vertices, triangles = reader.getMeshData()
        self.assertTrue((vertices == self.vertices).all())
        self.assertTrue((triangles == self.triangles).all())

    def test_3_read_group_data(self):
        """Read the group data from med file created in test 1.
        """
        reader = MEDTR3Reader(self.test_file)

        for group in self.groups:
            if group == "g3":
                # Group g3 is a collection of element ids, hence we expand it
                # to full mask
                file_group = reader.getGroupArray(group)
                _t = np.zeros(file_group.shape)
                _t[self.groups[group]] = 1
                self.assertTrue((file_group == _t).all())

            else:
                self.assertTrue((self.groups[group] == reader.getGroupArray(group)).all())

    def test_4_read_field_data(self):
        """Read the field data from med file created in test1.
        """
        reader = MEDTR3Reader(self.test_file)

        a = reader.getField('Test', 0)

        self.assertTrue(np.allclose(a, self.y))

        a = reader.getField('Test', 1)
        self.assertTrue(np.allclose(a, 2*self.y))

        a = reader.getField('Test', 2)
        self.assertTrue(np.allclose(a, 3*self.y))

        a = reader.getField('Test 2', 0)
        self.assertTrue(np.allclose(a, self.y2))

        a = reader.getField('Test 2', 1)
        self.assertTrue(np.allclose(a, self.y2*2))

    def test_4_Mesh(self):
        """Tests if the l2g.mesh.Mesh class operates well.
        """
        mesh = Mesh(self.test_file)

        vertices, triangles = mesh.getMeshData()
        self.assertTrue((vertices == self.vertices).all())
        self.assertTrue((triangles == self.triangles).all())

        for group in self.groups:
            if group == "g3":
                # Group g3 is a collection of element ids, hence we expand it
                # to full mask
                file_group = mesh.getGroupArray(group)
                _t = np.zeros(file_group.shape)
                _t[self.groups[group]] = 1
                self.assertTrue((file_group == _t).all())

            else:
                self.assertTrue((self.groups[group] == mesh.getGroupArray(group)).all())

        a = mesh.getField('Test', 0)

        self.assertTrue(np.allclose(a, self.y))

        a = mesh.getField('Test', 1)
        self.assertTrue(np.allclose(a, 2*self.y))

        a = mesh.getField('Test', 2)
        self.assertTrue(np.allclose(a, 3*self.y))

        a = mesh.getField('Test 2', 0)
        self.assertTrue(np.allclose(a, self.y2))

        a = mesh.getField('Test 2', 1)
        self.assertTrue(np.allclose(a, self.y2*2))

if __name__ == '__main__':
    unittest.main()