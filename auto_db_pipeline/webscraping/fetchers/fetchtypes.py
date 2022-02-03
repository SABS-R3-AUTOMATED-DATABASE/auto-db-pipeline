"""
Attempts to obtain the pmid and pmc for a paper.
"""
import warnings
import re
import requests
from metapub import PubMedFetcher
from metapub.exceptions import MetaPubError
from .fetchtext import get_html

class FetchTypes:
    """
    Object for fetching the types PMID and PMC.
    """

    def __init__(self, doi):
        self.doi = doi
        self.pmid, self.pmc = None, None
        self.found_metapub = None
        self.invalid_doi = False

    def __call__(self):
        if not self.doi:
            return
        self.try_metapub()
        self.try_pmc_redirect()
        self.try_manual_pubmed()

    def try_metapub(self):
        """
        Uses `metapub` library to attempt to obtain the PMID and PMC for a paper.
        """
        try:
            article_fetch = PubMedFetcher().article_by_doi(self.doi)
            self.found_metapub = True
            self.pmid = article_fetch.pmid
            self.pmc = article_fetch.pmc

        except MetaPubError:
            self.found_metapub = False

        except AttributeError:
            self.found_metapub = False
            self.invalid_doi = True

    def try_pmc_redirect(self):
        """
        Try to get the PMC by the redirect from the doi url.
        """
        if self.pmc:
            return
        request = requests.get(FetchTypes.get_url_pmc(self.doi))
        if request.status_code == 404:
            return
        pattern = re.compile(r"(PMC\d+)(/)")
        pattern_search = re.search(pattern, request.url)
        if pattern_search:
            self.pmc = pattern_search.group(1)

    def try_manual_pubmed(self):
        """
        Attempts to find a paper by manually querying the doi in pubmed.
        Uses patterns on the website to obtain the PMID and/or PMC.
        """
        url_search_pubmed = "https://pubmed.ncbi.nlm.nih.gov/?term=" + self.doi
        pubmed_html = get_html(url_search_pubmed)

        # Try to grab the PMID
        try:
            start_idx = pubmed_html.index("pmid:")
            # Length of "mid:" plus max length of PMID plus 1
            end_idx = start_idx + 4 + 8 + 1
            pmid_match = pubmed_html[start_idx:end_idx]
            if not self.pmid:
                self.pmid = re.search(r"\d+", pmid_match).group()
        except ValueError:  # substring not found
            start_idx = None  # useful for the next step when we check against start_idx

        if self.pmc:
            return

        # Try to grab the PMC ID
        try:
            self.pmc = re.findall(r"PMC\d{7}", pubmed_html)[0][3:]
            if start_idx:
                # Get index of the above and check that it is not far from pmid
                pmc_idx = pubmed_html.index("PMC" + self.pmc)
                if pmc_idx > start_idx+16:
                    setattr(self, "pmc_wrong_loc", True)
                    warnings.warn("PMC in wrong location.")

        except IndexError:  # list index out of range
            pass

    @staticmethod
    def get_url_pmc(doi):
        """Redirects to PMC"""
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/doi/{doi}"
