"""
Get proteininformatics from papers.
"""
from fetchtext import get_soup, get_text, get_html
from papertypes import Doi, Pmid, Pmc

class Paper:
    """
    Obtain the data on a paper.
    Input: a dictionary of a paper's metadata (`paper_data`)
    """

    def __init__(self, paper_data: dict):
        """
        Initialize and run the class to store the data in the instance of paper.
        """
        self.title = paper_data['title']
        self.idoi = paper_data['doi']
        self.authors = paper_data['authors']

        self.doi, self.pmid, self.pmc = Doi(), Pmid(), Pmc()

        self._clear_caches()


    def _clear_caches(self):
        """Clear the caches."""
        funcs_to_clear = (get_soup, get_text, get_html)
        for func in funcs_to_clear:
            func.cache_clear()
