import unittest
from Bio import Entrez

import mine_genbank


class TestMineGenbank(unittest.TestCase):
    '''class to test mine_genbank.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

        # initially search one entry
        self.search = 'QQX23546[Accession]'
        self.genbank = mine_genbank.GenbankSearch(self.search)
        self.genbank.get_number_of_entries()
        self.genbank.get_entries()
        # example of entires to test the functions
        self.entries_list = [{'GBSeq_references': [
                                 {'GBReference_authors': 'Author 1',
                                  'GBReference_title': 'Title 1'}],
                              'GBSeq_definition': 'HEAVY, SARS-cov, S3-F23'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 1',
                                  'GBReference_title': 'Title 1'}],
                              'GBSeq_definition': 'alpha, covid-19, RBD, A12'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 2',
                                  'GBReference_title': 'Title 2, merS-cOv'}],
                              'GBSeq_definition': 'kappa'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 2',
                                  'GBReference_title': 'Title 2, merS-cOv'}],
                              'GBSeq_definition': 'omicron'},
                             {'GBSeq_references': [
                                 {'GBReference_authors': 'Author 2',
                                  'GBReference_title': 'Title 1'}],
                              'GBSeq_definition': 'LAMBDA, ABD123'}]

        # list of ids of genbank entries, 1-2 & 3-4 are vh-vl pairings
        # 5-6 & 7-8 are from the same publication but not pairings
        self.ids = ['7WDF_F', '7WDF_E', '7WD7_f', '7WD7_e', '7WD7_C', '7WCZ_B',
                    'CAB4059600', 'CAB4065462']

        entires_handle = Entrez.efetch(db='protein', id=self.ids,
                                       rettype="gb", retmode="xml")
        self.entries_gb = Entrez.read(entires_handle)

    def test_get_number_of_entries(self):
        self.assertEqual(self.genbank.number_of_entries, 1)

    def test_get_entries(self):
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

    def test_classify_vh_vl(self):
        # overwrite entires with example entries
        self.genbank.entries = self.entries_list
        self.genbank.classify_vh_vl()
        # check if VHs and VLs in examples were classified correctly
        entries = self.genbank.entries
        self.assertEqual(entries[0]['chain'], 'heavy_chain')
        self.assertEqual(entries[1]['chain'], 'heavy_chain')
        self.assertEqual(entries[2]['chain'], 'light_chain')
        self.assertEqual(entries[3]['chain'], 'unassigned')
        self.assertEqual(entries[4]['chain'], 'light_chain')

    def test_find_antigen(self):
        # overwrite entires with example entries
        self.genbank.entries = self.entries_list
        self.genbank.find_antigen()
        # check if the antigen of the example entries was detected correctly
        entries = self.genbank.entries
        self.assertEqual(entries[0]['antigen'], 'SARS-CoV')
        self.assertEqual(entries[1]['antigen'],
                         'SARS-CoV-2, Spike protein RBD')
        self.assertEqual(entries[2]['antigen'], 'MERS-CoV')
        self.assertEqual(entries[3]['antigen'], 'MERS-CoV')
        self.assertEqual(entries[4]['antigen'], 'not determined')

    def test_find_fragment_id(self):
        # overwrite entires with example entries
        self.genbank.entries = self.entries_list
        self.genbank.find_fragment_id()
        # check if the fragment id of the example entries was determined
        # correctly
        entries = self.genbank.entries
        self.assertEqual(entries[0]['fragment_id'], 'S3-F23')
        self.assertEqual(entries[1]['fragment_id'], 'A12')
        self.assertEqual(entries[4]['fragment_id'], 'ABD123')
        with self.assertRaises(KeyError):
            entries[2]['fragment_id']

    def test_group_by_author_title(self):
        # overwrite entires with example entries
        self.genbank.entries = self.entries_list
        self.genbank.classify_vh_vl()
        self.genbank.group_by_author_title()
        grouped_entries = self.genbank.grouped_entries

        # check if the example entires were correctly grouped
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

        # check that all groups differ in either title or authors
        self.assertFalse(title_g1 == title_g2 and author_g1 == author_g2)
        self.assertFalse(title_g1 == title_g3 and author_g1 == author_g3)
        self.assertFalse(title_g2 == title_g3 and author_g2 == author_g3)

    def test_pair_vh_vl(self):
        # overwrite entires with genbank entries
        self.genbank.entries = self.entries_gb
        self.genbank.classify_vh_vl()
        self.genbank.find_antigen()
        self.genbank.find_fragment_id()
        self.genbank.group_by_author_title()
        self.genbank.pair_vh_vl()

        # check that the total number of entries after pairing equals
        # the total number of entires dowloaded from genbank
        total_entries = len(self.genbank.paired_entries) * 2
        # multiplied by two as each pairing has one heavy and light chain
        for group in self.genbank.unpaired_entries:
            total_entries += len(group)
        self.assertEqual(total_entries, len(self.ids))

        # check that the number of entires in each category are correct
        self.assertEqual(len(self.genbank.paired_entries), 2)
        self.assertEqual(len(self.genbank.unpaired_entries), 2)
        self.assertEqual(len(self.genbank.unpaired_entries[0]), 2)
        self.assertEqual(len(self.genbank.unpaired_entries[1]), 2)

        for entry in self.genbank.paired_entries:
            # check that the chains are correct in the dictionary
            self.assertEqual(entry['heavy_chain']['chain'], 'heavy_chain')
            self.assertEqual(entry['light_chain']['chain'], 'light_chain')

            # check that the pairings are correct
            if entry['heavy_chain']['GBSeq_locus'] == '7WDF_F':
                self.assertTrue(
                    entry['light_chain']['GBSeq_locus'] == '7WDF_E')
            if entry['heavy_chain']['GBSeq_locus'] == '7WD7_f':
                self.assertTrue(
                    entry['light_chain']['GBSeq_locus'] == '7WD7_e')


if __name__ == '__main__':
    unittest.main()
