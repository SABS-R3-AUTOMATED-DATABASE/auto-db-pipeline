"""
Object representations of the paper types: Doi, pmid, pmc.
"""
import requests
from .fetchers.fetchtext import get_text
from .proteinids.interface import Interface

class Doi:
    """
    Doi class.
    """
    def __init__(self, doi, authors, journal):
        if not doi:
            return
        self.doi = doi
        self.journal = journal
        self.set_url(journal)
        self.interface = Interface(doi, authors, self.paper_text)
        self.interface()

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
    """
    Pmid class.
    """
    def __init__(self, pmid, doi, authors):
        if not pmid:
            return
        self.pmid = pmid
        self.url = self.url_pmid
        self.interface = Interface(doi, authors, self.paper_text, pmid)
        self.interface()

    @property
    def url_pmid(self):
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"

    @property
    def paper_text(self):
        """Default paper text."""
        if self.url:
            return get_text(self.url)


class Pmc:
    """
    Pmc class.
    """
    def __init__(self, pmc, doi, authors):
        if not pmc:
            return
        self.pmc = pmc
        self.url = self.url_pmc
        self.interface = Interface(doi, authors, self.paper_text)
        self.interface()

    @property
    def url_pmc(self):
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{self.pmc}/"

    @property
    def paper_text(self):
        """Default paper text."""
        if self.url:
            return get_text(self.url)
