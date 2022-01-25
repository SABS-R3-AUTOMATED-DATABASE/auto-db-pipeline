"""
Unit tests of functions involved in the web-scraping.
"""
import unittest
import pdb_checker
import scrape_support


SKIP_SLOW_TESTS = True

class TestScraping(unittest.TestCase):
    """
    Unit tests of functions in the scrape support file.
    """

    def test_filter_pdb_id(self):
        """
        Unit test of filter_pdb function.
        """
        sample_paper_text = """(6VXX, 4IO4)
                            7LF4. PDB ID: 7THH!
                            "7SJ7-"
                            5ABC. 
                            .5abc, 
                            5126VXX5
                            """

        correct_pdb_tuple = ('6VXX', '4IO4', '7LF4', '7THH', '7SJ7', '5ABC', '5abc')

        filtered_pdb_tuple = tuple(scrape_support.filter_pdb_id(sample_paper_text))

        self.assertTupleEqual(correct_pdb_tuple, filtered_pdb_tuple)


    def test_filter_genbank_protein_id(self):
        """
        Unit test of filter_genbank function.
        """
        sample_paper_text = """
                            (e.g., aaa98665.1 will change to AAA98665.2). 
                            ABC12345f, (AAA98665.25). GenBank: 
                            QGA67030.1, ABC12345.2A, FABC12345.2A
                            """

        correct_genbank_id_tuple = ('aaa98665.1', 'AAA98665.2', 'AAA98665.25', 'QGA67030.1')
        filtered_genbank_id_tuple = tuple(scrape_support.filter_genbank_protein_id(sample_paper_text))

        self.assertTupleEqual(correct_genbank_id_tuple, filtered_genbank_id_tuple)


    def test_pdb_check(self):
        """
        Unit tests of pdb_check funciton.
        """

        sample_paper_texts_true = ["here is a (protein data bank-", "\nPDBID"]
        for sample_paper_text in sample_paper_texts_true:
            self.assertTrue(scrape_support.pdb_check(sample_paper_text))

        sample_paper_texts_true = ["protein data", "hello"]
        for sample_paper_text in sample_paper_texts_true:
            self.assertFalse(scrape_support.pdb_check(sample_paper_text))


    @unittest.skipIf(condition=SKIP_SLOW_TESTS, reason="This test has passed and is slow.")
    def test_get_top_authors(self):
        """
        Tests that the output of the get_top_authors function is correct for 
        PDBs: 6VXX.
        """

        pdb_id = "6VXX"
        top_num_authors = 3

        correct_top_authors = ('Walls', 'Park', 'Tortorici')

        checker = pdb_checker.PDBChecker()
        checker_top_authors = checker.get_top_authors(pdb_id, top_num_authors, verbose=False)

        for correct_author, checker_author in zip(correct_top_authors, checker_top_authors):
            self.assertEqual(correct_author.casefold(), checker_author.casefold())

if __name__ == "__main__":
    unittest.main()
