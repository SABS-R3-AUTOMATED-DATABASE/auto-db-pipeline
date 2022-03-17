"""
Representation of a paper.
Called by papers2ids.py
"""
from .fetchers.fetch_types import FetchTypes
from .fetchers.fetch_text import get_soup, get_text, get_html
from .paper_types import Doi, Pmid, Pmc

BIORXIV_NAME = 'biorxiv'
class Paper:
    """
    Obtain the data on a paper.
    Input: a dictionary of a paper's metadata (`paper_data`)"""

    def __init__(self, paper_data: dict):
        """
        Initialize and run the class to store the data in the instance of paper.
        """
        self.title = paper_data['title']
        self.journal = paper_data['journal']
        self.authors = paper_data['authors']
        self.doi = paper_data['doi']
        self.date = paper_data['date']
        self.source = paper_data['source']
        self.paper_types = {}

    def __repr__(self):
        """Representation."""
        return str(self.output)

    @property
    def output(self):
        """Get the output of the paper, for retrieval."""
        attrs = ('title', 'journal', 'authors', 'date', 'source')
        representation = {attr: getattr(self, attr) for attr in attrs}
        if not self.paper_types:
            return representation
        representation.update(self._output_paper_types)
        return representation

    def __enter__(self):
        """Enter the context manager and return the paper object."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context manager and clear the caches."""
        if exc_type:
            print(self.__dict__)
            print(exc_type)
            print(exc_value)
            print(traceback)
            print("\n")

        Paper._clear_caches()

    def __bool__(self):
        """If the doi does not exist, we can't get anything
        (with current methods)."""
        return bool(self.doi)

    def __call__(self):
        """Run the web-scraping and get the information."""
        print('calling paper')
        if not self:
            return
        self.fetch_paper_types()
        self.paper_type_interface()

    @property
    def _output_paper_types(self):
        """Get the output of the paper types, attach
        paper types that do not exist."""
        output_ = {}
        for type_name, paper_type in self.paper_types.items():
            output_[type_name] = paper_type.output
        return {'paper_types': output_}

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
    def _clear_caches():
        """Clear the caches."""
        funcs_to_clear = (get_soup, get_text, get_html)
        for func in funcs_to_clear:
            func.cache_clear()
