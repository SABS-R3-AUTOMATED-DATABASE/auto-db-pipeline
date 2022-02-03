"""
Get proteininformatics from papers.
"""
from .webscraping.fetchers.fetchtypes import FetchTypes
from .webscraping.fetchers.fetchtext import get_soup, get_text, get_html
from .webscraping.papertypes import Doi, Pmid, Pmc

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

        self.doi, self.pmid, self.pmc = Doi(), Pmid(), Pmc()
        self.fetch = FetchTypes(doi_string)
        self._initialize_types(doi_string)


    def __call__(self):
        """
        Run the web-scraping and get the information.
        """
        self.doi()
        self.pmid()
        self.pmc()
        self._clear_caches()

    def _initialize_types(self, doi_string):
        """
        Initialize the types by running the FetchTypes class and
        recording its results.
        """
        self.fetch()

        self.doi.set_doi(doi_string, self.journal)
        self.pmid.set_pmid(self.fetch.pmid)
        self.pmc.set_pmc(self.fetch.pmc)

    def _clear_caches(self):
        """
        Clear the caches.
        """
        funcs_to_clear = (get_soup, get_text, get_html)
        for func in funcs_to_clear:
            func.cache_clear()
