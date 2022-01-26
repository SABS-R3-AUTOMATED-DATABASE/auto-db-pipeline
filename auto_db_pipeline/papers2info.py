"""
File docstring.
"""
import warnings
import re
import requests
from metapub import PubMedFetcher
from pdb_checker import PDBChecker
from support import filter_pdb_id, filter_genbank_protein_id
from support import get_soup, get_html, get_text, pdb_check, genbank_check

class Paper:
    """
    Obtain the data on a paper.

    Input: a dictionary of a paper's metadata (`paper_data`) which comes
    from a query using the 'paperscraper' library. See the below:
    https://github.com/PhosphorylatedRabbits/paperscraper

    Interacts with a `checker`, which is an instance of the `PDBChecker` class.

    * Attempts to find the final correct re-directed url of the paper
      from the doi (to get the full text)
    * Attempts to obtain the PMID and PMC for a paper and their urls (for the full text).
    * Attempts to scrape the PDB IDs of a paper from the full text urls (web-pages).
    * Note that PDB ID web-scraping occurs only if 'PDB' is in the paper text. See comment # remove
    """

    def __init__(self, paper_data: dict, verbose=True):
        """
        Initialize and run the class to store the data in the instance of paper.
        """
        self.title, self.authors, self.date, self.abstract, self.journal, self.doi = paper_data.values()
        self.verbose = verbose
        self.pmid, self.pmc = None, None

        self.pdb_checker = PDBChecker()
        self.set_web_types()
        self.retrieve_paper()
        self.get_pdb()
        self.get_genbank_protein()
        self.clear_caches()


    def clear_caches(self):
        """
        Clear the caches.
        """
        funcs_to_clear = (get_soup, get_text, get_html)
        for func in funcs_to_clear:
            func.cache_clear()


    def set_web_types(self):
        """
        Set the attributes for the web-types, e.g. PMID, PMC.
        """
        self.web_type_names = ['doi', 'pmid', 'pmc']
        self.url_type_names = [f'url_{web_type}' for web_type in self.web_type_names]

        for web_type_name, url_type_name in zip(self.web_type_names, self.url_type_names):
            setattr(self, url_type_name, None)

        self.pdb_in_paper, self.genbank_in_paper, self.actual_pdbs, self.poss_gb_p_ids = {}, {}, {}, {}
        for web_type_name in self.web_type_names:
            self.pdb_in_paper[web_type_name] = False
            self.genbank_in_paper[web_type_name] = False
            self.actual_pdbs[web_type_name] = None
            self.poss_gb_p_ids[web_type_name] = None



    def loop_texts(self):
        """
        Loop through the available webpage texts for scraping functions.
        """
        for web_type_name, url_type_name in zip(self.web_type_names, self.url_type_names):
            url_type = getattr(self, url_type_name)
            if url_type:
                paper_text = get_text(url_type)
                yield web_type_name, paper_text


    def get_genbank_protein(self):
        """
        Retrieve the possible GenBank Protein IDs mentioned in each paper text version.
        """
        if self.verbose:
            print("Retrieving GenBank Protein IDs and checking mentions in papers...")
        for web_type_name, paper_text in self.loop_texts():
            self.genbank_in_paper[web_type_name] = genbank_check(paper_text)
            self.poss_gb_p_ids[web_type_name] = filter_genbank_protein_id(paper_text)
        if self.verbose:
            print(self.poss_gb_p_ids, end="\n\n")


    def get_pdb(self):
        """
        Retrieve the actual PDB IDs mentioned in each paper text version.
        These can be combined later.
        Also checks if 'pdb' or 'protein data bank' is mentioned in the text of the paper.
        (Not case sensitive.)
        """
        if self.verbose:
            print("Retrieving PDB IDs and checking mentions in papers...")
        for web_type_name, paper_text in self.loop_texts():
            self.pdb_in_paper[web_type_name] = pdb_check(paper_text)
            poss_pdbs = filter_pdb_id(paper_text)
            actual_pdbs = self.pdb_checker.get_actual(poss_pdbs)
            self.actual_pdbs[web_type_name] = actual_pdbs
        if self.verbose:
            print(self.actual_pdbs, end="\n\n")
            print(self.pdb_in_paper, end="\n\n")


    @property
    def any_pdb_in_paper(self):
        """
        Returns a boolean of whether 'PDB' or 'protein data bank' was found in any of the
        paper texts.
        """
        return any(self.pdb_in_paper.values())


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

        if self.verbose:
            print("DOI exists. Attempting to retrieve...", end='\n\n')

        if self.verbose:
            print("Getting url from the doi...", end="\n\n")
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


    def get_url_doi(self):
        """
        Get the url that results from the doi.org approach.
        """
        setattr(self, "url_doi", requests.get('https://doi.org/' + self.doi, allow_redirects=True).url)

        if ("medRxiv" in self.journal) or ("bioRxiv" in self.journal):
            suffix = ".full"
        else:  # room to add more alterations
            suffix = ""
        self.url_doi += suffix


    def get_url_pmid(self):
        # Add "/" to check for redirects later
        setattr(self, "url_pmid", "https://pubmed.ncbi.nlm.nih.gov/" + self.pmid + "/")


    def get_url_pmc(self):
        # Add "/" to check for redirects later
        setattr(self, "url_pmc", "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC" + self.pmc + "/")
