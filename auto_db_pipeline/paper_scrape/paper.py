"""
Representation of a ppaper.
Called by papers2ids.py
"""
from .fetchers.fetch_types import FetchTypes
from .fetchers.fetch_text import get_soup, get_text, get_html
from .paper_types import Doi, Pmid, Pmc
from ..protein_scrape.pdb_scrape import PdbID
from ..protein_scrape.genbank_scrape import GenBankID

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
        self.journal = Paper.set_journal(source, paper_data)
        self.authors = paper_data['authors']
        self.doi = paper_data['doi']
        self.date = paper_data['date']
        self.source = source
        self.paper_types = {}

    def __repr__(self):
        to_keep = ('title', 'journal', 'authors', 'date', 'source')
        representation = {attr: self.__dict__[attr] for attr in to_keep}
        if not self.paper_types:
            return representation
        representation.update(self._repr_paper_types)
        return representation

    def __enter__(self):
        """Enter the context manager and return the paper object."""
        return self

    def __exit__(self):
        """Exit the context manager and clear the caches."""
        Paper._clear_caches()

    def __bool__(self):
        """If the doi does not exist, we can't get anything
        (with current methods)."""
        return bool(self.doi)

    def __call__(self):
        """Run the web-scraping and get the information."""
        if not self:
            return
        self.fetch_paper_types()
        self.paper_type_interface()

    @property
    def _repr_paper_types(self):
        """Get the representation of the paper types, attach
        paper types that do not exist."""
        repr_ = {}
        for type_name, paper_type in self.paper_types.items():
            repr_[type_name] = repr(paper_type)
        return {'paper_types': repr_}

    def fetch_paper_types(self):
        """Fetch the pmid and pmc identifiers using the FetchTypes class.
        Then initialize the paper type objects."""
        fetch = FetchTypes(self.doi)
        self.paper_types['doi'] = Doi(self.doi, self.authors, self.journal)
        self.paper_types['pmid'] = Pmid(fetch.pmid, self.doi, self.authors)
        self.paper_types['pmc'] = Pmc(fetch.pmc, self.doi, self.authors)

    def paper_type_interface(self):
        for paper_type in self.paper_types.values():
            if paper_type:
                paper_type.interface()

    @staticmethod
    def set_journal(source, paper_data):
        """Set the journal, which is bioRxiv when the source is bioRxiv."""
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
