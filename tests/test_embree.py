import unittest

from l2g.comp.core import PyEmbreeAccell
import numpy as np

class TestEmbree(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        """Use class method setUpClass so that the embreeObject is only
        instantiated once. Using the setUp function creates a new instance
        of Embree object every time.
        """
        cls.embreeObj = PyEmbreeAccell()

    def test_1_add_first_mesh(self):
        """Add a mesh or geometry to Embree object. It's geomId should be
        0.
        """
        vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float32)
        triangles = np.array([1,2,3], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 0)

    def test_2_add_second_mesh(self):
        """Adds a second mesh or geometry to Embree object. It's geomId should
        be incremented by 1.
        """
        vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float32)
        vertices += 1.0
        triangles = np.array([1,2,3], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 1)

    def test_3_remove_last_mesh(self):
        """Removing first commited geometry.
        """
        self.assertTrue(self.embreeObj.deleteMesh(0))

    def test_4_add_mesh_after_removal(self):
        """After removing the first geometry with geometry ID = 0, we add
        another mesh to see if the geomId is 0.
        """
        vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float32)
        vertices += 1.0
        triangles = np.array([1,2,3], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 0)

if __name__ == '__main__':
    unittest.main()