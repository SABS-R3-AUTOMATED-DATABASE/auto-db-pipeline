import unittest
from Bio import Entrez

import mine_genbank


class TestMineGenbank(unittest.TestCase):
    '''class to test mine_genbank.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

        self.search = 'QQX23546[Accession]'
        self.genbank = mine_genbank.GenbankSearch(self.search)
        self.genbank.get_number_of_entries()
        self.genbank.get_entries()

    def test_get_number_of_entries(self):
        self.assertEqual(self.genbank.number_of_entries, 1)

    def test_get_entries(self):
        self.assertIsInstance(self.genbank.entries, list)
        self.assertIsInstance(self.genbank.entries[0], dict)
        self.assertEqual(self.genbank.entries[0]['GBSeq_moltype'], 'AA')

        self.sequence = 'evqllesgggvvqpgrslrlscaasgftfssygmhwvrqapgkglewvav' +\
                        'iwydgsnkyyadsvkgrftisrdnskntlylqmnslraedtavyycardg' +\
                        'sgsyywgsldywgqgtlvtvss'
        self.assertEqual(self.genbank.entries[0]['GBSeq_sequence'],
                         self.sequence)


if __name__ == '__main__':
    unittest.main()
