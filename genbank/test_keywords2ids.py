import unittest
from unittest.mock import MagicMock
from unittest.mock import patch, mock_open
from Bio import Entrez
import json

import keywords2ids


class TestKeywords2ids(unittest.TestCase):
    '''class to test keywords2ids.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

        self.search = 'QQX23546[Accession]'
        self.genbank = keywords2ids.GenbankSearch(self.search)
        self.genbank.get_number_of_entries()
        self.genbank.get_ids()

    def test_get_number_of_entries(self):
        self.assertEqual(self.genbank.number_of_entries, 1)

    def test_get_ids(self):
        # check that the formatting of the dowloaded ids is correct
        self.assertIsInstance(self.genbank.ids, list)
        # check that only one id was retrieved
        self.assertEqual(len(self.genbank.ids), 1)
        # check id is in a string format
        self.assertIsInstance(self.genbank.ids[0], str)
        # check the id is correct
        self.assertEqual(self.genbank.ids[0], '1960820478')

    def test_save_to_json(self):
        json.dump = MagicMock()
        m = mock_open(read_data='version= 1.0.0')
        with patch('builtins.open', m):
            self.genbank.save_to_json()
        m.assert_called_once_with('genbank/data/id_list.json', 'w')
        json.dump.assert_called_once()


if __name__ == '__main__':
    unittest.main()
