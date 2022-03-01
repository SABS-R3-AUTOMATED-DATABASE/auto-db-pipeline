import unittest
from unittest.mock import Mock, MagicMock
import json
from pandas import DataFrame
import pandas as pd

import info2csv


class TestInfo2csv(unittest.TestCase):
    '''class to test info2csv.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)

        self.paired = [{"heavy_chain": {"GBSeq_create-date": "23-FEB-2022",
                                        "GBSeq_definition": "Chain H",
                                        "GBSeq_source": "Macaca mulatta",
                                        "GBSeq_sequence": "cbacba",
                                        "GBSeq_references":
                                            [{"GBReference_authors":
                                                ["He,W.", "Name, A."],
                                              "GBReference_title":
                                                "Neutralizing antibodies",
                                              "GBReference_journal":
                                                "Unpublished"}],
                                        "chain": "H",
                                        "antigen": "SARS-CoV",
                                        "fragment_id": "K398.22"},
                        "light_chain": {"GBSeq_create-date": "23-FEB-2022",
                                        "GBSeq_definition": "Chain L",
                                        "GBSeq_accession-version": "7TP4_L",
                                        "GBSeq_source": "Macaca mulatta",
                                        "GBSeq_sequence": "abcabc",
                                        "GBSeq_references":
                                        [{"GBReference_authors":
                                            ["He,W.", "Name, A."],
                                          "GBReference_title":
                                              "Neutralizing antibodies",
                                          "GBReference_journal":
                                              "Unpublished"}],
                                        "chain": "L",
                                        "antigen": "SARS-CoV",
                                        "fragment_id": "K398.22"}}]

        self.unpaired = [[{"GBSeq_create-date": "23-FEB-2022",
                           "GBSeq_definition": "Chain H",
                           "GBSeq_source": "Macaca mulatta",
                           "GBSeq_sequence": "cbacba",
                           "GBSeq_references":
                               [{"GBReference_authors":
                                   ["He,W.", "Name, A."],
                                 "GBReference_title":
                                   "Neutralizing antibodies",
                                 "GBReference_journal":
                                   "Unpublished"}],
                           "chain": "H",
                           "antigen": "SARS-CoV",
                           "fragment_id": "K398.22"}]]

        # initialise with above examples for paired and unpaires seqs
        json.load = Mock(side_effect=[self.paired, self.unpaired])
        self.genbank = info2csv.PopulateDatabase()

    def test_create_df(self):
        df = self.genbank.create_df()
        # test type dataframe
        self.assertIsInstance(df, DataFrame)
        # test dataframe is empty
        self.assertEqual(len(df), 0)

    def test_get_refs(self):
        ref = self.genbank.get_ref(self.paired[0]['heavy_chain'])
        # test if reference is correct
        true_ref = 'He,W. et al.; Neutralizing antibodies; Unpublished; 2022.'
        self.assertEqual(ref, true_ref)

    def test_prepare_paired_row(self):
        row_dict = self.genbank.prepare_paired_row(self.paired[0])
        self.assertEqual(row_dict['Name'], 'K398.22')
        self.assertEqual(row_dict['Origin'], 'Macaca mulatta')
        self.assertEqual(row_dict['Binds_to'], 'SARS-CoV')

    def test_prepare_unpaired_row(self):
        row_dict = self.genbank.prepare_unpaired_row(self.unpaired[0][0])
        self.assertEqual(row_dict['Name'], 'K398.22')
        self.assertEqual(row_dict['Origin'], 'Macaca mulatta')
        self.assertEqual(row_dict['Binds_to'], 'SARS-CoV')

    def test_populate_db_paried(self):
        self.genbank.populate_db_paired()
        # test type dataframe
        self.assertIsInstance(self.genbank.df_p, DataFrame)
        # test dataframe has 1 entry
        self.assertEqual(len(self.genbank.df_p), 1)
        # test that entry is correct
        self.assertEqual(self.genbank.df_p.iloc[0, 0], 'K398.22')

    def test_populate_db_unparied(self):
        self.genbank.populate_db_unpaired()
        # test type dataframe
        self.assertIsInstance(self.genbank.df_up, DataFrame)
        # test dataframe has 1 entry
        self.assertEqual(len(self.genbank.df_up), 1)
        # test that entry is correct
        self.assertEqual(self.genbank.df_up.iloc[0, 0], 'K398.22')

    def test_filter_duplicates(self):
        # duplicate the paired sequences
        self.genbank.paired_seqs.append(self.paired[0])
        self.genbank.populate_db_paired()
        self.genbank.filter_duplicates_paried()

        # test type dataframe
        self.assertIsInstance(self.genbank.df_p, DataFrame)
        # test dataframe has 1 entry, second one removed
        self.assertEqual(len(self.genbank.df_p), 1)
        # test that entry is correct
        self.assertEqual(self.genbank.df_p.iloc[0, 0], 'K398.22')

    def test_save_to_csv(self):

        self.genbank.populate_db_paired()
        pd.DataFrame.to_csv = MagicMock()
        path = 'path'
        self.genbank.save_to_csv(self.genbank.df_p, path)
        pd.DataFrame.to_csv.assert_called_once_with(path)


if __name__ == '__main__':
    unittest.main()
