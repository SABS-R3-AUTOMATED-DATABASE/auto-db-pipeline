"""
This file contains the code to support the webscraping, such as in scrape.py.
One example is to retrieve potential/possible PDB IDs from a url
(e.g. the url of a paper).
and return those potential PDB IDs as a list.
"""
import functools
import re
import requests
from bs4 import BeautifulSoup

PAPER_SOUP_CACHE_TIMEOUT = 10

@functools.lru_cache(maxsize=PAPER_SOUP_CACHE_TIMEOUT)
def get_soup(url):
    """
    Gets all the text from the html of a url of a paper.
    Returns the prettified text as a string.

    Requires the following libraries:
    * requests
    * BeautifulSoup from bs4
    """
    # Get the response to the request to the url
    request = requests.get(url, allow_redirects=True)
    # Get the beautiful soup
    soup = BeautifulSoup(request.content, 'html.parser')
    return soup


@functools.lru_cache(maxsize=PAPER_SOUP_CACHE_TIMEOUT)
def get_text(url: str) -> str:
    """
    Convert the html in the beautiful soup into text (type string)
    """
    soup = get_soup(url)
    return soup.get_text()


@functools.lru_cache(maxsize=PAPER_SOUP_CACHE_TIMEOUT)
def get_html(url: str) -> str:
    """
    Return a prettified version of the html.
    """
    soup = get_soup(url)
    return soup.prettify()



def pdb_check(paper_text):
    """
    Checks if 'pdb' or 'protein data bank' is mentioned in a paper text.
    There are no English words that have PDB so we
    can feel confident that it is only finding references to PDB IDs.
    """
    pattern = re.compile(r'(pdb|protein data bank)', re.IGNORECASE)
    pdb_in_paper = bool(pattern.search(paper_text))
    return pdb_in_paper


def genbank_check(paper_text):
    """
    Checks if 'pdb' or 'protein data bank' is mentioned in a paper text.
    There are no English words that have PDB so we
    can feel confident that it is only finding references to PDB IDs.
    """
    pattern = re.compile(r"""(genbank|National Genetic Sequence Data Base)""",
                        re.IGNORECASE)
    genbank_in_paper = bool(pattern.search(paper_text))
    return genbank_in_paper


def standardize(filter_id_func):
    """
    Apply the standardization to a list of protein IDs.
    """
    def wrapper(*args, **kwargs):
        return standardize_ids(filter_id_func(*args, **kwargs))

    return wrapper


@standardize
def filter_pdb_id(paper_text: str) -> list:
    """
    Takes a paper text and returns a list of possible PDBs, i.e,  alphanumeric sequences
    that follow the rules of PDB IDs, given below:

    1. All characters are alphabetical or numeric
    1. Length of four.
    2. First charater is numeric.
    3. First character is in the range of 1-9 inclusive (greater than 0).

    These rules are taken from [protopedia](https://proteopedia.org/wiki/index.php/PDB_code).
    """
    pattern = re.compile(r"""
                        [^A-Z0-9]  # no alphanumeric character before
                        (  # begin group of interest
                        [1-9]  # first character of PDB is 1-9
                        [A-Z0-9]{3}  # remaining 3 characters of PDB are alphanumeric
                        )  # end group of interest
                        [^A-Z0-9]  # no alphanumeric character after
                        """, flags=re.IGNORECASE | re.VERBOSE)  # not case sensitive
    pdb_list = pattern.findall(paper_text)
    return pdb_list


@standardize
def filter_genbank_protein_id(paper_text: str) -> list:
    """
    Takes a paper_text (webpage text) and returns a list of potential
    GenBank protein IDs (alphanumeric paper_texts that follow the GenBank
    ID rules, explained below).

    The rules for protein GenBank IDs can be found here:
    * https://www.ncbi.nlm.nih.gov/genbank/sequenceids/
      "The protein IDs contain three letters followed by five digits,
       a period, and a version number."
    * https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
      "Protein IDs consist of three letters followed by five digits,
      a dot, and a version number."

    """
    pattern = re.compile(r"""
                        [^A-Z0-9]  # no alphanumeric character before
                        (  # begin group of interest
                        [A-Z]{3}  # three letters
                        \d{5}  # five digits
                        \.  # one period
                        \d+  # 1 or more digits for the version number
                        )  # end group of interest
                        [^A-Z]  # no alphabetic character after
                        """, flags=re.IGNORECASE | re.VERBOSE)  # not case sensitive
    genbank_protein_id_list = pattern.findall(paper_text)
    return genbank_protein_id_list


def standardize_ids(id_list: list) -> list:
    """
    Takes a list of possible protein IDs and does the following:
    1. Converts them to upper case.
    2. Removes duplicates.
    3. Sorts in alphabetical-numerical order.

    Misc note:
    We can remove any of these steps and this function will not throw an error.
    """
    # 1. Convert to upper case
    id_list = list(map(lambda id_: id_.upper(), id_list))
    # 2. Remove duplicates
    id_list = list(set(id_list))
    # 3. Sort alphabetically / numerically
    id_list.sort()
    return id_list

