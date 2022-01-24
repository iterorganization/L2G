import unittest

from l2g.comp.core import PyEmbreeAccell, PyFLT
import numpy as np

"""Check
"""

class TestModule(unittest.TestCase):

    def test_flt_object_creation(self):
        with self.assertRaises(Exception):
            try:
                obj = PyFLT()
            except:
                pass
            else:
                raise Exception

    def test_embree_object_creation(self):
        with self.assertRaises(Exception):
            try:
                obj = PyEmbreeAccell()
            except:
                pass
            else:
                raise Exception

    def test_assign_embree_to_flt(self):
        # We do not need to populate the embree object with data, the upper
        # function does that.
        with self.assertRaises(Exception):
            try:
                embree_obj = PyEmbreeAccell()
                flt_obj = PyFLT()
                flt_obj.applyRT(embree_obj)
            except:
                pass
            else:
                raise Exception

    def test_embree_populate_with_float_data(self):
        """Create 2 meshes and commit it to EmbreeObject. What is tested here
        are the actual geometry IDs that Embree outputs. Per Embree
        documentation the ID starts with 0 and then incremented for each
        geometry.
        """
        obj = PyEmbreeAccell()

        vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float32)
        triangles = np.array([1,2,3], dtype=np.uint32)
        geomId = obj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 0)

    def test_embree_populate_with_double_data(self):
        obj = PyEmbreeAccell()

        vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float64)
        triangles = np.array([1,2,3], dtype=np.uint32)
        geomId = obj.commitMesh(vertices, triangles)

        self.assertEqual(geomId, 0)

        # Now add another mesh
        vertices = np.array(vertices, copy=True)
        vertices += 1
        triangles = np.array(triangles, copy=True)
        geomId = obj.commitMesh(vertices, triangles)
        self.assertEqual(geomId, 1)

        # Now delete the last mesh and add a new mesh.

    def test_embree_populate_with_1D_python_list(self):
        """Pass points and triangles as 1D python lists.
        """
        obj = PyEmbreeAccell()

        vertices = [0, 0, 0, 1, 0, 0, 0, 1, 0]
        triangles = [1, 2, 3]

        # We just see if there is an error.
        with self.assertRaises(Exception):
            try:
                geomId = obj.commitMesh(vertices, triangles)
            except Exception as e:
                # If an exception is raised don't do anything.
                pass
            else:
                raise Exception

    def test_embree_populate_with_2D_python_list(self):
        """Pass points and triangles as 2D python lists.

        In this case we test is we have a list of points and a list of
        triangles, where every element is a list of 3 values.
        """
        obj = PyEmbreeAccell()

        vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        triangles = [[1, 2, 3]]

        # We just see if there is an error.
        with self.assertRaises(Exception):
            try:
                geomId = obj.commitMesh(vertices, triangles)
            except Exception as e:
                # If an exception is raised don't do anything.
                pass
            else:
                raise Exception

    def test_embree_populate_with_incorrect_type(self):
        """Pass points from different type than what Embree expects. I.e.,
        Embree expects float for vertices and unsigned integers for indexes
        of the triangles. Cython function should correctly recast the data
        to the appropriate data type.
        """
        obj = PyEmbreeAccell()

        vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float64)
        triangles = np.array([1,2,3], dtype=np.int64)

        # We just see if there is an error.
        with self.assertRaises(Exception):
            try:
                geomId = obj.commitMesh(vertices, triangles)
            except Exception as e:
                # If an exception is raised don't do anything.
                pass
            else:
                raise Exception


if __name__ == '__main__':
    unittest.main()
