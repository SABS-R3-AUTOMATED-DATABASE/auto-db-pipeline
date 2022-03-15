import unittest
from unittest.mock import MagicMock
from unittest.mock import patch, mock_open, Mock
import json

import json_combination


class TestIds2ids(unittest.TestCase):
    '''class to test ids2ids.py'''

    def __init__(self, *argv, **kwarg):
        super().__init__(*argv, **kwarg)

        # intitiallise with test genbank id
        self.id_test = ['A', 'B', 'C']
        self.id_papers_test = ['C', 'D', 'E']
        json.load = Mock(side_effect=[self.id_test, self.id_papers_test])
        m = mock_open(read_data='version= 1.0.0')
        with patch('builtins.open', m):
            self.genbank = json_combination.Combination()

    def test_combine_lists(self):
        self.genbank.combine_lists()
        self.assertEqual(len(self.genbank.ids), 5)

    def test_save_to_json(self):
        json.dump = MagicMock()
        m = mock_open(read_data='version= 1.0.0')
        with patch('builtins.open', m):
            self.genbank.save_to_json()
        m.assert_called_once_with('genbank/data/id_combined_list.json', 'w')
        json.dump.assert_called_once()


if __name__ == '__main__':
    unittest.main()
