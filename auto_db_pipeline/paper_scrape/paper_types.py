"""
Representations of the paper types: Doi, pmid, pmc.
Called by paper.py

Note that we can optionally choose to include the paper text
in the representations.
"""
import requests
from .fetchers.fetch_text import get_text
from .types_interface import Interface

class Doi:
    """Doi class."""
    def __init__(self, doi, authors, journal):
        self.doi = doi
        if not self:
            return
        print('getting doi')
        self.journal = journal
        self.set_url(journal)
        self.interface = Interface(doi, authors, self.paper_text)

    def __bool__(self):
        """Return whether the paper type can be considered to exist."""
        return bool(self.doi)

    def __repr__(self):
        """Representation."""
        return str(self.output)

    @property
    def output(self):
        """Output of the doi paper type for retrieval."""
        representation = {'doi': self.doi}
        if not self:
            return representation
        representation.update(self.interface.output)
        return representation

    def set_url(self, journal):
        """Set the url for the Doi class."""
        if "medRxiv" in journal:
            self.url = self.url_medrxiv
        elif "bioRxiv" in journal:
            self.url = self.url_biorxiv
        elif journal in ('Viruses', 'Vaccines'):
            url_redirect = Doi.get_redirect(self.url_doi)
            self.url = url_redirect + "/htm"
        else:
            self.url = self.url_doi

    @property
    def url_biorxiv(self):
        return f"https://www.biorxiv.org/content/{self.doi}.full"

    @property
    def url_medrxiv(self):
        return f"https://www.medrxiv.org/content/{self.doi}.full"

    @property
    def url_pmc(self):
        """Redirects to PMC"""
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/doi/{self.doi}"

    @property
    def url_doi(self):
        return f"https://doi.org/{self.doi}"

    @staticmethod
    def get_redirect(url):
        request = requests.get(url, allow_redirects=True)
        return request.url

    @property
    def paper_text(self):
        """Default paper text."""
        if self.url:
            return get_text(self.url)

class Pmid:
    """Pmid class."""
    def __init__(self, pmid, doi, authors):
        self.pmid = pmid
        if not self:
            return
        print('getting pmid')
        self.url = self.url_pmid
        self.interface = Interface(doi, authors, self.paper_text, pmid)

    def __bool__(self):
        """Return whether the paper type can be considered to exist."""
        return bool(self.pmid)

    def __repr__(self):
        """Representation."""
        return str(self.output)

    @property
    def output(self):
        """Output of the pmid paper type for retrieval."""
        representation = {'pmid': self.pmid}
        if not self:
            print('reached here')
            return None
        representation.update(self.interface.output)
        return representation

    @property
    def url_pmid(self):
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"

    @property
    def paper_text(self):
        """Default paper text."""
        if self.url:
            return get_text(self.url)

class Pmc:
    """Pmc class."""
    def __init__(self, pmc, doi, authors):
        self.pmc = pmc
        if not self:
            return
        print('getting pmc')
        self.url = self.url_pmc
        self.interface = Interface(doi, authors, self.paper_text)

    def __bool__(self):
        """Return whether the paper type can be considered to exist."""
        return bool(self.pmc)

    def __repr__(self):
        """Representation."""
        return str(self.output)

    @property
    def output(self):
        """Output of the pmc paper type for retrieval."""
        representation = {'pmc': self.pmc}
        if not self:
            return None
        representation.update(self.interface.output)
        return representation

    @property
    def url_pmc(self):
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{self.pmc}/"

    @property
    def paper_text(self):
        """Default paper text."""
        if self.url:
            return get_text(self.url)

class PaperType:
    """General paper type class."""

    def __init__(self, type_id):
        self.type_id = type_id
        if not self:
            return

    def __bool__(self):
        """Return whether the paper type can be considered to exist."""
        return bool(self.type_id)

    # Can probably remove doi and pmid from types_interface
