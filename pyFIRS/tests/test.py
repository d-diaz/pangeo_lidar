import os
import unittest
import numpy as np
from pyFIRS.utils import listlike

this_dir = os.path.dirname(__file__)


class TestUtils(unittest.TestCase):

    def test_listlike(self):
        """Checks whether listlike function returns appropriate boolean
        responses."""

        self.assertTrue(listlike([1,2]))  # list
        self.assertTrue(listlike((1,2)))  # tuple
        self.assertTrue(listlike(np.array((1,2.5,4))))  # array
        self.assertTrue({'a':1})  # dictionary
        self.assertFalse(listlike(1))  # single number
        self.assertFalse(listlike('string'))  # string


if __name__ == '__main__':
    unittest.main()
