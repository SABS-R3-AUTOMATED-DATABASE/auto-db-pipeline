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
        self.entries_list = [{'GBSeq_references': [
                                 {'GBReference_authors': 'Author 1',
                                  'GBReference_title': 'Title 1'}],
                              'GBSeq_definition': 'heavy'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 1',
                                  'GBReference_title': 'Title 1'}],
                              'GBSeq_definition': 'alpha'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 2',
                                  'GBReference_title': 'Title 2'}],
                              'GBSeq_definition': 'kappa'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 2',
                                  'GBReference_title': 'Title 2'}],
                              'GBSeq_definition': 'omicron'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 2',
                                  'GBReference_title': 'Title 1'}],
                              'GBSeq_definition': 'lambda'}]

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

    def test_classify_vh_vl(self):
        self.genbank.entries = self.entries_list
        self.genbank.classify_vh_vl()
        entries = self.genbank.entries
        self.assertEqual(entries[0]['chain'], 'heavy_chain')
        self.assertEqual(entries[1]['chain'], 'heavy_chain')
        self.assertEqual(entries[2]['chain'], 'light_chain')
        self.assertEqual(entries[3]['chain'], 'unassigned')
        self.assertEqual(entries[4]['chain'], 'light_chain')

    def test_group_by_author_title(self):
        self.genbank.entries = self.entries_list
        self.genbank.classify_vh_vl()
        self.genbank.group_by_author_title()
        grouped_entries = self.genbank.grouped_entries

        self.assertEqual(len(grouped_entries), 3)
        self.assertEqual(len(grouped_entries[0]), 2)

        g1 = grouped_entries[0][0]['GBSeq_references'][0]
        g2 = grouped_entries[1][0]['GBSeq_references'][0]
        g3 = grouped_entries[2][0]['GBSeq_references'][0]
        author_g1 = g1['GBReference_authors']
        author_g2 = g2['GBReference_authors']
        author_g3 = g3['GBReference_authors']
        title_g1 = g1['GBReference_title']
        title_g2 = g2['GBReference_title']
        title_g3 = g3['GBReference_title']

        self.assertFalse(title_g1 == title_g2 and author_g1 == author_g2)
        self.assertFalse(title_g1 == title_g3 and author_g1 == author_g3)
        self.assertFalse(title_g2 == title_g3 and author_g2 == author_g3)


if __name__ == '__main__':
    unittest.main()
