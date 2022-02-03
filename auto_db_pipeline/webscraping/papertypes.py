"""
Object representations of the paper types: Doi, pmid, pmc.
"""
import requests
from .fetchers.fetchtext import get_text
from .proteinids.proteinids import exists_mention_of_id_type, get_instances_of_id
from .proteinids.proteinids import id_checking, id_finding
from .proteinids.pdbinterface import PdbInterface

class PaperType:
    """
    General class for a type of paper text, (e.g. doi, pmid, and pmc).
    """
    def __init__(self):
        self.mention_ids = {id_name: None for id_name in id_checking}
        self.possible_ids = {id_name: None for id_name in id_finding}
        self.actual_ids = {id_name: None for id_name in id_finding}
        self.url = None  # Default url

    def __bool__(self):
        """Fancy way to check if the instance of this PaperType exists."""
        return bool(self.url)

    @property
    def paper_text(self):
        """Default paper text."""
        if self.url:
            return get_text(self.url)

    def __call__(self):
        self.get_all_id_mention()
        self.get_all_id_possible()
        self.get_pdb_id_actual()

    def get_all_id_mention(self):
        for id_name in id_checking:
            self._get_mention(id_name)

    def get_all_id_possible(self):
        for id_name in id_finding:
            self._get_possible(id_name)

    def get_pdb_id_actual(self):
        self.actual_ids['pdb_id'] = PdbInterface().get_actual(self.possible_ids['pdb_id'])

    def _get_mention(self, id_name):
        regex = id_checking[id_name]
        self.mention_ids[id_name] = exists_mention_of_id_type(self.paper_text, regex)

    def _get_possible(self, id_name):
        regex = id_finding[id_name]
        self.possible_ids[id_name] = get_instances_of_id(self.paper_text, regex)


class Doi(PaperType):
    """
    Doi class.
    """
    def __init__(self):
        super().__init__()
        self.doi = None

    def set_doi(self, doi, journal):
        if not doi:
            return
        self.doi = doi
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


class Pmid(PaperType):
    """
    Pmid class.
    """
    def __init__(self):
        super().__init__()
        self.pmid = None

    def set_pmid(self, pmid):
        if not pmid:
            return
        self.pmid = pmid
        self.url = self.url_pmid

    @property
    def url_pmid(self):
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"


class Pmc(PaperType):
    """
    Pmc class.
    """
    def __init__(self):
        super().__init__()
        self.pmc = None

    def set_pmc(self, pmc):
        if not pmc:
            return
        self.pmc = pmc
        self.url = self.url_pmc

    @property
    def url_pmc(self):
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{self.pmc}/"
