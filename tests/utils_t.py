import sys
sys.path.insert(1,"./")

from typing import List
import unittest
import numpy as np

from src import Utils as u

class UtilsTest(unittest.TestCase):

    def setUp(self):
        self.data = u.read_dw_data("data/white_dwarf_data.csv")
        
    def test_read_dw_data_first_element(self):
        self.assertEqual(np.linalg.norm(self.data[0][:] - np.array([7.455, 0.301])), 0)
    
    def test_read_dw_data_last_element(self):
        self.assertEqual(np.linalg.norm(self.data[-1][:] - np.array([8.139, 0.665])), 0)

    def test_read_dw_data_data_type_ndarray(self):
        self.assertIsInstance(self.data, np.ndarray)

    def test_read_dw_data_data_type_list(self):
        self.data = u.read_dw_data("data/white_dwarf_data.csv", out_type="list")
        self.assertIsInstance(self.data, List)

    def test_read_dw_data_shape_0(self):
        self.assertEqual(self.data.shape[0], 378)

    def test_read_dw_data_shape_1(self):
        self.assertEqual(self.data.shape[1], 2)
if __name__ == "__main__":

    unittest.main()
