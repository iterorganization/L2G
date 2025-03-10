import unittest

from l2g.external.embree import PyEmbreeAccell
import numpy as np

"""Test the Embree API interface of loading/removing meshes and testing if
ray intersection works.
"""

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
        vertices = np.array([0,0,0,
                             1,0,0,
                             0,1,0], dtype=np.float32)
        triangles = np.array([0,1,2], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 0)

    def test_2_add_second_mesh(self):
        """Adds a second mesh or geometry to Embree object. It's geomId should
        be incremented by 1.
        """
        vertices = np.array([1,0,0,
                             2,0,0,
                             1,1,0], dtype=np.float32)
        triangles = np.array([0,1,2], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 1)

    def test_3_remove_first_mesh(self):
        """Removing first commited geometry.
        """
        self.assertTrue(self.embreeObj.deleteMesh(0))

    def test_4_add_mesh_after_removal(self):
        """After removing the first geometry with geometry ID = 0, we add
        another (the same) mesh to see if the geomId is 0.
        """
        vertices = np.array([0,0,0,
                             1,0,0,
                             0,1,0], dtype=np.float32)
        triangles = np.array([0,1,2], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 0)

    def test_5_add_mesh_with_name(self):
        """Add a mesh with name 'test' ask the Embree object if it contains a
        mesh with that name.
        """
        vertices = np.array([2,0,0,
                             3,0,0,
                             2,1,0], dtype=np.float32)
        triangles = np.array([0,1,2], dtype=np.uint32)
        geomId = self.embreeObj.commitMesh(vertices, triangles, 'test')

        isIn = self.embreeObj.isMeshWithNameIn('test')
        self.assertTrue(isIn)

    def test_6_test_intersection(self):
        """Test intersection by putting the ray perpendicular on the loaded
        triangle.
        """

        self.embreeObj.castInfRay(0.33, 0.33, -1, 0, 0, 1)
        self.assertTrue(self.embreeObj.checkIfHit())
        self.assertEqual(self.embreeObj.returnGeomId(), 0)
        self.assertEqual(self.embreeObj.returnPrimId(), 0)

        self.embreeObj.castInfRay(0.33, 0.33, -1, 0, 0, -1)
        self.assertFalse(self.embreeObj.checkIfHit())

        # Cast it into the third triangle or triangle named test.
        self.embreeObj.castInfRay(2.33, 0.33, -1, 0, 0, 1)
        self.assertTrue(self.embreeObj.checkIfHit())
        self.assertEqual(self.embreeObj.returnGeomId(), 2)

if __name__ == '__main__':
    unittest.main()
