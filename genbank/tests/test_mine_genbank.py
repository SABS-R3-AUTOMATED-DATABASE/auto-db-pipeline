import unittest
# import numpy as np
# import numpy.testing as npt
from unittest.mock import patch
# import pytest
import pandas as pd
from Bio import Entrez

from genbank import mine_genbank


class TestMineGenbank(unittest.TestCase):
    '''class to test mine_genbank.py'''

    def __init__(self):
        self.search = 'QQX23546[Accession]'
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

    def test_initialize(self):

        with patch.object(pd, 'read_csv') as mock_read_csv:
            CoV = Keywords.Keywords('./test/path.csv')
            mock_read_csv.assert_called_with('./test/path.csv')
            print(CoV)


if __name__ == '__main__':
    unittest.main()
