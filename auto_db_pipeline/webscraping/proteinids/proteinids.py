"""
This file contains the code to support the webscraping, such as in scrape.py.
One example is to retrieve potential/possible PDB IDs from a url
(e.g. the url of a paper).
and return those potential PDB IDs as a list.
"""
import re

id_checking = {'genbank': r"(genbank|National Genetic Sequence Data Base)",
                'pdb': r"(pdb|protein data bank)",
                'accession': r"(accession)",
                'protein': r"(protein|antibody|antibodies)",
                'nucleotide': r"nucleotide",
                'geninfo': r"(geninfo|gi number)",
                'refseq': r"(refseq|reference sequence)"
            }

id_finding = {
    'pdb_id': r"""
                        \b  # starts at word boundary
                        (
                            [1-9]  # first character of PDB is 1-9
                            [A-Z0-9]{3}  # remaining 3 characters of PDB are alphanumeric
                        )
                        [^A-Z0-9]  # no alphanumeric character after, but possibly _ for chain
                        """,
    'genbank_protein_id': r"""
                        \b  # starts at word boundary
                        (
                            [A-Z]{3}  # three letters
                            \d{5}  # five digits
                            \.  # one period
                            \d+  # one or more digits for the version number
                        )
                        [^A-Z]  # no alphabetic character after
                        """,
    'genbank_protein_accession': r"""
                        \b  # begins at word boundary
                        (
                            [A-Z]{3}  # three letters
                            (?:  # non-caputuring group
                                \d{5}|\d{7}  # five or seven digits
                            )
                        )
                        \b  # ends at word boundary
                        """,
    'genbank_nucleotide_accession': r"""
                        \b # begins at word boundary
                        (
                            (?:  # one letter and five digits
                                [A-Z]{1}\d{5}
                            )
                            |  # or
                            (?:
                                [A-Z]{2}  # two letters and
                                (?: \d{6}|\d{8})  # six or eight digits
                            )
                        )
                        \b  # ends at word boundary
                        """,
    'refseq_id': r"""
                        \b # begins at word boundary
                        [A-Z]{2}  # two letters
                        _  # underscore
                        \d{6,}  # six or more digits
                        \b # ends at word boundary
                        """,
    'geninfo_identifier': r"""
                        \b  # begins at word boundary
                        \d+  # at least one digit
                        \b  # ends at word boundary
                        """}

def _standardize(filter_id_func):
    """
    Apply the standardization to a list of protein IDs.
    """
    def wrapper(*args, **kwargs):
        return _standardize_ids(filter_id_func(*args, **kwargs))
    return wrapper

@_standardize
def get_instances_of_id(paper_text: str, regex: str) -> list:
    """
    Gets all the mentions of an ID in a paper, using a regex
    string in `id_finding`.
    """
    # Not case sensitive
    pattern = re.compile(regex, flags=re.IGNORECASE | re.VERBOSE)
    instances = pattern.findall(paper_text)
    return instances

def exists_mention_of_id_type(paper_text: str, regex: str) -> bool:
    """
    Checks whether an ID type is mentioned in a paper, using a regex
    string in `id_checking`.
    """
    pattern = re.compile(regex, re.IGNORECASE)
    id_in_paper = bool(pattern.search(paper_text))
    return id_in_paper

def _standardize_ids(id_list: list) -> list:
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
