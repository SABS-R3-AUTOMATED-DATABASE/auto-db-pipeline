"""
Attempts to obtain the Doi, pmid, and pmc for a paper.
"""
import warnings
import re
import requests
from metapub import PubMedFetcher
from metapub.exceptions import MetaPubError
from fetchtext import get_html
class FetchTypes:

    def __init__(self, doi):
        self.doi = doi
        self.pmid, self.pmc = None, None
        self.found_metapub = None
        self.invalid_doi = False

    def get_pmid_and_pmc(self, doi):
        self.try_metapub()
        if self.found_metapub:
            return

    def try_metapub(self):
        """
        Uses `metapub` library to attempt to obtain the PMID and PMC for a paper.
        """
        try:
            article_fetch = PubMedFetcher().article_by_doi(self.doi)
            if self.pmc and self.pmid:
                self.found_metapub = True
            self.pmid = article_fetch.pmid
            self.pmc = article_fetch.pmc

        except MetaPubError:
            self.found_metapub = False

        except AttributeError:
            self.invalid_doi = True


    def try_pmc_redirect(self):
        request = requests.get(FetchTypes.get_url_pmc(self.doi))
        if request.status_code == 404:
            return
        pattern = re.compile(r"(PMC\d+)(/)")
        pattern_search = re.search(pattern, request.url)
        if pattern_search:
            self.pmc = pattern_search.group(1)


    @staticmethod
    def get_url_pmc(doi):
        """Redirects to PMC"""
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/doi/{doi}"



    def try_pubmed(self):
        """
        Attempts to find a paper by manually querying the doi in pubmed.
        Uses patterns on the website to obtain the PMID and/or PMC.
        """

        url_search_pubmed = "https://pubmed.ncbi.nlm.nih.gov/?term=" + self.doi

        pubmed_html = get_html(url_search_pubmed)

        try:
            start_idx = pubmed_html.index("pmid:")
            # Length of "mid:" plus max length of PMID plus 1
            end_idx = start_idx + 4 + 8 + 1
            pmid_match = pubmed_html[start_idx:end_idx]
            setattr(self, "pmid", pmid_match.split(':')[1].split(',')[0])
        except ValueError:  # substring not found
            if self.verbose:
                print("PMID not found via pubmed webscrape.")
            setattr(self, "pmid", None)
            start_idx = None  # useful for the next step when we check against start_idx

        try:
            setattr(self, "pmc", re.findall(r"PMC\d{7}", pubmed_html)[0][3:])
            if start_idx:
                # Get index of the above and check that it is not far from pmid
                pmc_idx = pubmed_html.index("PMC" + self.pmc)
                if pmc_idx > start_idx+16:
                    warnings.warn("PMC in wrong location.")

        except IndexError:  # list index out of range
            if self.verbose:
                print("PMC not found via pubmed webscrape.")
            setattr(self, "pmc", None)

        if (not self.pmid) and not (self.pmc):
            setattr(self, "found_pubmed", False)
        else:
            setattr(self, "found_pubmed", True)

        if not getattr(self, "found_pubmed"):
            if self.verbose:
                print("Pubmed webscrape failed. No PMID or PMC available.", end='\n\n')
        else:
            if self.verbose:
                print("Pubmed webscrape succceeded.", end='\n\n')

