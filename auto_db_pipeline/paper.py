"""
Get proteininformatics from papers.
"""
from webscraping.fetchers.fetchtypes import FetchTypes
from webscraping.fetchers.fetchtext import get_soup, get_text, get_html
from webscraping.papertypes import Doi, Pmid, Pmc
from webscraping.proteinids.idtypes import PdbID, GenBankID

BIORXIV_NAME = 'biorxiv'

class Paper:
    """
    Obtain the data on a paper.
    Input: a dictionary of a paper's metadata (`paper_data`)"""

    def __init__(self, paper_data: dict, source: str):
        """
        Initialize and run the class to store the data in the instance of paper.
        """
        self.title = paper_data['title']
        self.authors = paper_data['authors']
        self.journal = Paper.set_journal(source, paper_data)
        self.doi = paper_data['doi']
        self.source = source
        self.paper_types = {}

    def __call__(self):
        """Run the web-scraping and get the information."""
        self.fetch_paper_types()
        self.interface()
        Paper._clear_caches()

    def fetch_paper_types(self):
        fetch = FetchTypes(self.doi)
        self.paper_types['doi'] = Doi(self.doi, self.authors, self.journal)
        self.paper_types['pmid'] = Pmid(fetch.pmid, self.doi, self.authors)
        self.paper_types['pmc'] = Pmc(fetch.pmc, self.doi, self.authors)

    def interface(self):
        for paper_type in self.paper_types.values():
            if paper_type:
                paper_type.interface()

    @staticmethod
    def set_journal(source, paper_data):
        if source == BIORXIV_NAME:
            return 'bioRxiv'
        else:
            return paper_data['journal']

    @staticmethod
    def _clear_caches():
        """Clear the caches."""
        funcs_to_clear = (get_soup, get_text, get_html,
                            GenBankID.get_handle, PdbID.get_handle)
        for func in funcs_to_clear:
            func.cache_clear()
