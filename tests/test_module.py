from tabnanny import verbose
import unittest

from l2g.external.tlas import PyTLAS
from l2g.external.flt import PyFLT
from l2g.external.equilibrium_analysis import EQA
from l2g.external.bicubic import PyBicubic
from l2g.external.bfgs_2d import PyBfgs2d
from l2g.external.rkf45 import PyRKF45FLT
import numpy as np


class TestModule(unittest.TestCase):

    def test_flt_object_creation(self):
        for _ in range(10):
            obj = PyFLT()

    def test_tlas_object_creation(self):
        for _ in range(10):
            obj = PyTLAS()

    def test_assign_tlas_to_flt(self):
        for _ in range(10):
            tlas_obj = PyTLAS()
            flt_obj = PyFLT()
            flt_obj.applyRT(tlas_obj)


    def test_PyBfgs2d_object_creation(self):
        for _ in range(10):
            obj = PyBfgs2d()
            del obj

    def test_PyRKF45FLT_object_creation(self):
        for _ in range(10):
            obj = PyRKF45FLT()
            del obj

    def test_EQA_object_creation(self):
        for _ in range(10):
            obj = EQA()
            del obj

if __name__ == '__main__':
    unittest.main()
