import unittest
# import numpy as np
# import numpy.testing as npt
from unittest.mock import patch
# import pytest
import pandas as pd

from keywords import Keywords


class TestKeywords(unittest.TestCase):
    '''class to test Keywords'''

    def test_initialize(self):

        with patch.object(pd, 'read_csv') as mock_read_csv:
            CoV = Keywords.Keywords('./test/path.csv')
            mock_read_csv.assert_called_with('./test/path.csv')
            print(CoV)


if __name__ == '__main__':
    unittest.main()
