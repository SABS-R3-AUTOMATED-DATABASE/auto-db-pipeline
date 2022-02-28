import unittest
from Bio import Entrez
import numpy

import evaluate_genbank_search


class TestMineGenbank(unittest.TestCase):
    '''class to test evaluate_genbank_searches.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

        self.keywords = 'unittest'
        self.protein_entries = [{
            'GBSeq_sequence': 'QVQLVQSGAEVKKPGESLKISCKGSGYSFTSYWIGWVRQMPGKG' +
                              'LEWMGIIYPGDSDTRYSPSFQGQVTISADKSISTAYLQWSSLKA' +
                              'SDTAMYYCARVGSYQPSFDYWGQGTLVTVSS'
            },
            {
            'GBSeq_sequence': 'ABC',
            },
            {
            'GBSeq_sequence': 'EIVMTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAP' +
                              'RLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYC' +
                              'QQRSNWPWTFGQGTKVEIK'
            }
        ]

        self.evaluation = evaluate_genbank_search.EvaluateGenbankSearch(
                          'genbank/data/CoV-AbDab_181021.csv',
                          self.protein_entries, self.keywords)
        self.evaluation(print_metrics=False, save_metrics=False)

    def test_prepare_variables(self):
        self.assertEqual(len(self.evaluation.aa_seqs), 3)
        self.assertIsInstance(self.evaluation.aa_seqs, numpy.ndarray)
        self.assertEqual(len(self.evaluation.VH_arr),
                         len(self.evaluation.VH_found))

    def test_produce_metrics(self):
        self.assertEqual(self.evaluation.total_seqs, 3)
        self.assertEqual(self.evaluation.sequences_in_covabdab, 2)
        self.assertEqual(self.evaluation.sequences_not_in_covabdab, 1)
        self.assertEqual(self.evaluation.match_rate, 2/3)
        self.assertEqual(len(self.evaluation.CovAbDab_stats),
                         self.evaluation.no_seqs_in_covabdab)
        self.assertEqual(self.evaluation.covabdab_vh_found, 1)
        self.assertEqual(self.evaluation.covabdab_vl_found, 1)
        self.assertEqual(self.evaluation.VH_VL_pairings_found, 0)
        self.assertEqual(self.evaluation.VH_or_VL_found, 2)
        self.assertAlmostEqual(self.evaluation.percentage_VH_or_VL_found,
                               2 * 100 / self.evaluation.no_seqs_in_covabdab)


if __name__ == '__main__':
    unittest.main()
