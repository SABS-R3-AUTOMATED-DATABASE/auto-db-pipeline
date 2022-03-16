import unittest
from unittest.mock import MagicMock
from unittest.mock import patch, mock_open
from Bio import Entrez
import json

import ids2protein


class TestIds2protein(unittest.TestCase):
    '''class to test ids2protein.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

        # intitiallise with test genbank id
        self.test_id = ['1960820478']
        json.load = MagicMock(return_value=self.test_id)
        self.genbank = ids2protein.ProteinRetrieval()

    def test_get_entries(self):
        self.genbank.get_entries()

        # check that the formatting of the dowloaded entry is correct
        self.assertIsInstance(self.genbank.entries, list)
        self.assertIsInstance(self.genbank.entries[0], dict)
        # check that the correct (protein) database was searched
        self.assertEqual(self.genbank.entries[0]['GBSeq_moltype'], 'AA')

        # check that the downloaded sequence is correct
        self.sequence = 'evqllesgggvvqpgrslrlscaasgftfssygmhwvrqapgkglewvav' +\
                        'iwydgsnkyyadsvkgrftisrdnskntlylqmnslraedtavyycardg' +\
                        'sgsyywgsldywgqgtlvtvss'
        self.assertEqual(self.genbank.entries[0]['GBSeq_sequence'],
                         self.sequence)

    def test_remove_junk(self):
        self.genbank.get_entries()
        self.genbank.remove_junk()

        # remove junk retains nine key-value pairs
        self.assertEqual(len(self.genbank.entries[0]), 9)

    def test_save_to_json(self):
        self.genbank.get_entries()
        self.genbank.remove_junk()

        json.dump = MagicMock()
        m = mock_open(read_data='version= 1.0.0')
        with patch('builtins.open', m):
            self.genbank.save_to_json()
        m.assert_called_once_with('genbank/data/protein_handles.json', 'w')
        json.dump.assert_called_once()


if __name__ == '__main__':
    unittest.main()
