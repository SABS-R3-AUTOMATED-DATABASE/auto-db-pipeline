"""
Unit tests of functions involved in the web-scraping.
"""
import unittest
from auto_db_pipeline.webscraping.proteinids.extractids import exists_mention, get_instances
from auto_db_pipeline.webscraping.proteinids.extractids import id_finding, id_checking

class TestScraping(unittest.TestCase):
    """
    Unit tests of functions in the scrape support file.
    """
    def test_get_instances_pdb(self):
        """
        Unit test of filter_pdb function.
        """
        sample_paper_text = """(6VXX_H, 4IO4)
                            7LF4. PDB ID: 7THH!
                            "7SJ7-"
                            5ABC.
                            .5abc,
                            5126VXX5
                            """
        correct_pdb_tuple = ('4IO4', '5ABC', '6VXX', '7LF4', '7SJ7', '7THH')
        regex = id_finding['pdb_id']
        filtered_pdb_tuple = tuple(get_instances(sample_paper_text, regex))
        self.assertTupleEqual(correct_pdb_tuple, filtered_pdb_tuple)

    def test_get_instances_genbank(self):
        """
        Unit test of filter_genbank function.
        """
        sample_paper_text = """(e.g., aaa98665.1 will change to AAA98665.2).
                            ABC12345f, (AAA98665.25). GenBank:
                            QGA67030.1, FABC13345.2
                            QGA67030
                            QGA67031F
                            """
        correct_genbank_id_tuple = ('AAA98665.1', 'AAA98665.2', 'AAA98665.25', 'QGA67030.1')
        regex = id_finding['genbank_protein_id']
        filtered_genbank_id_tuple = tuple(get_instances(sample_paper_text, regex))
        self.assertTupleEqual(correct_genbank_id_tuple, filtered_genbank_id_tuple)

    def test_exists_mention(self):
        """
        Unit tests of pdb_check funciton.
        """
        sample_paper_texts_true = ["here is a (protein data bank-", "\nPDBID"]
        regex = id_checking['pdb']
        for sample_paper_text in sample_paper_texts_true:
            self.assertTrue(exists_mention(sample_paper_text, regex))

        sample_paper_texts_true = ["protein data", "hello"]
        for sample_paper_text in sample_paper_texts_true:
            self.assertFalse(exists_mention(sample_paper_text, regex))

if __name__ == "__main__":
    unittest.main()
