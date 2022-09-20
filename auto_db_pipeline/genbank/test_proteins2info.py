import unittest
from unittest.mock import MagicMock
from unittest.mock import patch, mock_open
import json
from Bio import Entrez

import proteins2info


class TestProteins2info(unittest.TestCase):
    '''class to test proteins2info.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

        self.test_entries = [{'GBSeq_references': [
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

        self.known_antigens = {'COVID': ['coronavirus', 'sars-cov', 'sars'],
                               'MERS': ['mers-cov', 'mers'],
                               'SARS-COV-2': ['sars-cov-2', 'covid-19'],
                               'SARS-COV-1': ['sars-cov-1'],
                               'Spike protein': ['spike', 'spike protein'],
                               'RBD': ['receptor binding domain', 'rbd']}

        json.load = MagicMock(return_value=self.test_entries)
        self.genbank = proteins2info.InfoRetrieval()

        # list of ids of genbank entries, 1-2 & 3-4 are vh-vl pairings
        # 5-6 & 7-8 are from the same publication but not pairings
        self.test_ids = ['7WDF_F', '7WDF_E', '7WD7_f', '7WD7_e', '7WD7_C',
                         '7WCZ_B', 'CAB4059600', 'CAB4065462']

        entires_handle = Entrez.efetch(db='protein', id=self.test_ids,
                                       rettype="gb", retmode="xml")
        self.entries_gb = Entrez.read(entires_handle)

        self.entries_nb = [{'GBSeq_definition': 'sybody'},
                           {'GBSeq_definition': 'antibody'},
                           {'GBSeq_definition': 'nanobody'}]

    def test_filter_AB_entries(self):
        # overwrite with test entries
        self.genbank.entries = self.entries_gb
        # filter keeps 4 entries: '7WDF_F', '7WDF_E', '7WD7_f', '7WD7_e'
        self.genbank.filter_AB_entries()
        # check if 4 entries were retatined
        self.assertEqual(len(self.genbank.entries), 4)
        # check if the correct entires are retained
        ids = []
        correct_ids = ['7WDF_F', '7WDF_E', '7WD7_f', '7WD7_e']
        for entry in self.genbank.entries:
            ids.append(entry['GBSeq_locus'])
        self.assertEqual(ids, correct_ids)

    def test_classify_vh_vl(self):
        self.genbank.classify_vh_vl()
        # check if VHs and VLs in examples were classified correctly
        entries = self.genbank.entries
        self.assertEqual(entries[0]['chain'], 'H')
        self.assertEqual(entries[1]['chain'], 'H')
        self.assertEqual(entries[2]['chain'], 'L')
        self.assertEqual(entries[3]['chain'], 'unassigned')
        self.assertEqual(entries[4]['chain'], 'L')

    def test_classify_vh_vl_anarci(self):
        # overwrite with first 4 test entries which are ABs
        # chains are H, L, L, H in this order
        self.genbank.entries = self.entries_gb[:4]
        self.genbank.classify_vh_vl_anarci()
        # check if VHs and VLs in examples were classified correctly
        entries = self.genbank.entries
        self.assertEqual(entries[0]['chain'], 'H')
        self.assertEqual(entries[1]['chain'], 'L')
        self.assertEqual(entries[2]['chain'], 'L')
        self.assertEqual(entries[3]['chain'], 'H')

    def test_find_antigen(self):
        self.genbank.find_antigen(self.known_antigens)
        # check if the antigen of the example entries was detected correctly
        entries = self.genbank.entries
        self.assertEqual(entries[0]['antigen'], ['COVID'])
        self.assertEqual(entries[1]['antigen'],
                         ['SARS-COV-2', 'RBD'])
        self.assertEqual(entries[2]['antigen'], ['MERS'])
        self.assertEqual(entries[3]['antigen'], ['MERS'])
        self.assertEqual(entries[4]['antigen'], 'unknown')

    def test_find_fragment_id(self):
        self.genbank.find_fragment_id()
        # check if the fragment id of the entries is correct
        entries = self.genbank.entries
        self.assertEqual(entries[0]['fragment_id'], 'S3-F23')
        self.assertEqual(entries[1]['fragment_id'], 'A12')
        self.assertEqual(entries[4]['fragment_id'], 'ABD123')
        with self.assertRaises(KeyError):
            entries[2]['fragment_id']

    def test_find_nanobodies(self):
        # use self.entries_nb to test nanobody functions
        self.genbank.entries = self.entries_nb
        self.genbank.find_nanobodies()

        # test correct classification
        entries = self.genbank.entries
        self.assertEqual(entries[0]['antibody_type'], 'nanobody')
        self.assertEqual(entries[1]['antibody_type'], 'antibody')
        self.assertEqual(entries[2]['antibody_type'], 'nanobody')

    def test_seperate_nanobodies(self):
        # use self.entries_nb to test nanobody functions
        self.genbank.entries = self.entries_nb
        self.genbank.find_nanobodies()
        self.genbank.seperate_nanobodies()

        # test correct classification
        self.assertEqual(len(self.genbank.entries), 1)
        self.assertEqual(len(self.genbank.nanobodies), 2)

    def test_group_by_publication(self):
        self.genbank.classify_vh_vl()
        self.genbank.group_by_publication()
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
        # overwrite entires with example entries
        self.genbank.entries = self.entries_gb
        self.genbank.classify_vh_vl()
        self.genbank.find_antigen(self.known_antigens)
        self.genbank.find_fragment_id()
        self.genbank.group_by_publication()
        self.genbank.pair_vh_vl()

        # check that the total number of entries after pairing equals
        # the total number of entires dowloaded from genbank
        total_entries = len(self.genbank.paired_entries) * 2
        for group in self.genbank.unpaired_entries:
            total_entries += len(group)
        self.assertEqual(total_entries, len(self.test_ids))

        # check that the number of entires in each category are correct
        self.assertEqual(len(self.genbank.paired_entries), 2)
        self.assertEqual(len(self.genbank.unpaired_entries), 2)
        self.assertEqual(len(self.genbank.unpaired_entries[0]), 2)
        self.assertEqual(len(self.genbank.unpaired_entries[1]), 2)

        for entry in self.genbank.paired_entries:
            # check that the chains are correct in the dictionary
            self.assertEqual(entry['heavy_chain']['chain'], 'H')
            self.assertEqual(entry['light_chain']['chain'], 'L')

            # check that the pairings are correct
            if entry['heavy_chain']['GBSeq_locus'] == '7WDF_F':
                self.assertTrue(
                    entry['light_chain']['GBSeq_locus'] == '7WDF_E')
            if entry['heavy_chain']['GBSeq_locus'] == '7WD7_f':
                self.assertTrue(
                    entry['light_chain']['GBSeq_locus'] == '7WD7_e')

    def test_get_chains_from_sabdab(self):
        entry = {'GBSeq_locus': '7RR0_B',
                 'GBSeq_sequence': 'qvqlvqsgpevkkpgtsvkvsckasgftfsssavqwvrqa' +
                                   'rgqrlewigwivvgsgnanyaprfqervtitrdmstntay' +
                                   'melsslrsedtavyycaapncsrtlcydgfnmwgqgtmvtv',
                 'chain': 'H'}

        result = self.genbank.get_chains_from_sabdab(entry)

        HC = 'QVQLVQSGPEVKKPGTSVKVSCKASGFTFSSSAVQWVRQARGQRLEWIGWIVVGSGNANYA' +\
             'PRFQERVTITRDMSTNTAYMELSSLRSEDTAVYYCAAPNCSRTLCYDGFNMWGQGTMVTV'
        LC = 'EIVLTQSPGSLSLSPGERATLSCRASQSVRSSYLGWYQQKPGQAPRLLIYGASSRATGIPD' +\
             'RFSGSGSETDFTLTISRLEPEDFAVYYCQQYDSSPWTFGQGTKVEI'
        self.assertEqual(HC, result[0])
        self.assertEqual(LC, result[1])

    def test_get_pairable_sequences_from_sabdab(self):
        self.genbank.get_chains_from_sabdab = MagicMock(return_value=(1, 2))
        # example of nested list input used by function
        unpaired_test = [[{'1': 1}, {'1': 1}], [{'1': 1}]]
        self.genbank.unpaired_entries = unpaired_test
        self.genbank.get_pairable_sequences_from_sabdab()
        # with example input Mock should be called 4 times
        self.assertEqual(self.genbank.get_chains_from_sabdab.call_count, 3)

    def test_pair_vh_vl_pdb(self):
        pdb_pariables_test = [{'chain': 'H',
                               'pdb_lc': 'ABC',
                               'GBSeq_sequence': 'dcba'},
                              {'chain': 'L',
                               'pdb_hc': 'CBA',
                               'GBSeq_sequence': 'abcd'}]

        self.genbank.paired_entries = []
        self.genbank.pdb_pairables = pdb_pariables_test
        self.genbank.pair_vh_vl_pdb()

        self.assertEqual(len(self.genbank.paired_entries), 1)
        test = self.genbank.paired_entries[0]['heavy_chain']['GBSeq_sequence']
        self.assertEqual(test, 'dcba')

    def test_save_to_json(self):
        # create onjects for paired, unpaired and nanobody entires
        self.genbank.paired_entries = []
        self.genbank.unpaired_entries = []
        self.genbank.nanobodies = []

        json.dump = MagicMock()
        m = mock_open(read_data='version= 1.0.0')
        with patch('builtins.open', m):
            self.genbank.save_to_json()
        self.assertEqual(m.call_count, 3)
        self.assertEqual(json.dump.call_count, 3)


if __name__ == '__main__':
    unittest.main()
