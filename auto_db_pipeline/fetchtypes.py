"""
Attempts to obtain the Doi, pmid, and pmc for a paper.
"""
import re
import warnings
from metapub import PubMedFetcher
from fetchtext import get_html


    def retrieve_paper(self):
        """
        Retrieves the paper full-text dois for webscraping.
        """
        if not self.doi:
            # Cannot proceed of no doi
            if self.verbose:
                print("DOI does not exist. Catastrophic failure.", end='\n\n')
            self.url_doi = None
            return  # empty return

        self.get_url_doi()

        if self.verbose:
            print("Trying metapub fetch...")
        self.try_metapub()

        if not getattr(self, "found_metapub"):
            # We only try the pubmed search `try_pubmed` if metapub at least partly fails.
            # Therefore we are assuming that if metapub succeeds, it will always
            # find the PMID and PMC if they exist. We should #test this.
            if self.verbose:
                print("Trying pubmed webscrape...")
            self.try_pubmed()

        if self.pmid:
            self.get_url_pmid()
        if self.pmc:
            self.get_url_pmc()

        if self.verbose:
            print("DOI:", self.doi)
            print("PMID:", self.pmid)
            print("PMC:", self.pmc)


    def try_metapub(self):
        """
        Uses `metapub` library to attempt to obtain the PMID and PMC for a paper.
        """
        try:
            article_fetch = PubMedFetcher().article_by_doi(self.doi)
            setattr(self, "found_metapub", True)
            self.pmid = article_fetch.pmid
            self.pmc = article_fetch.pmc

            if self.verbose:
                if self.pmid and self.pmc:
                    print("Metapub fetch completely succeeded. PMID and PMC found.", end='\n\n')
        except Exception:
            setattr(self, "found_metapub", False)
            if self.verbose:
                print("Metapub fetch at least partly failed.", end='\n\n')


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

