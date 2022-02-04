"""
Get proteininformatics from papers.
"""
from .webscraping.fetchers.fetchtypes import FetchTypes
from .webscraping.fetchers.fetchtext import get_soup, get_text, get_html
from .webscraping.papertypes import Doi, Pmid, Pmc

# pylint: disable=too-few-public-methods

class Paper:
    """
    Obtain the data on a paper.
    Input: a dictionary of a paper's metadata (`paper_data`)
    """

    def __init__(self, paper_data: dict, journal=None):
        """
        Initialize and run the class to store the data in the instance of paper.
        """
        self.title = paper_data['title']
        self.authors = paper_data['authors']
        if journal:
            self.journal = journal
        else:
            self.journal = paper_data['journal']

        doi_string = paper_data['doi']

        self.fetch = FetchTypes(doi_string)

        self.doi = Doi(doi_string, self.authors, self.journal)
        self.pmid = Pmid(self.fetch.pmid, doi_string, self.authors)
        self.pmc = Pmc(self.fetch.pmc, doi_string, self.authors)


    def __call__(self):
        """
        Run the web-scraping and get the information.
        """
        self.fetch()
        self.doi.interface()
        self.pmid.interface()
        self.pmc.interface()

        Paper._clear_caches()


    @staticmethod
    def _clear_caches():
        """
        Clear the caches.
        """
        funcs_to_clear = (get_soup, get_text, get_html)
        for func in funcs_to_clear:
            func.cache_clear()
