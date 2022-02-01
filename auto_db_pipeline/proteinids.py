"""
This file contains the code to support the webscraping, such as in scrape.py.
One example is to retrieve potential/possible PDB IDs from a url
(e.g. the url of a paper).
and return those potential PDB IDs as a list.
"""
import re

id_checking = {'genbank': r"""(genbank|National Genetic Sequence Data Base)""",
                'pdb': r"""(pdb|protein data bank)"""}

id_finding = {'genbank_protein_id': r"""
                                    [^A-Z0-9]  # no alphanumeric character before
                                    (
                                    [1-9]  # first character of PDB is 1-9
                                    [A-Z0-9]{3}  # remaining 3 characters of PDB are alphanumeric
                                    )
                                    [^A-Z0-9]  # no alphanumeric character after
                                    """,
                'pdb_id': r"""
                            [^A-Z0-9]  # no alphanumeric character before
                            (
                            [A-Z]{3}  # three letters
                            \d{5}  # five digits
                            \.  # one period
                            \d+  # 1 or more digits for the version number
                            )
                            [^A-Z]  # no alphabetic character after
                            """
                                }


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


def _standardize(filter_id_func):
    """
    Apply the standardization to a list of protein IDs.
    """
    def wrapper(*args, **kwargs):
        return _standardize_ids(filter_id_func(*args, **kwargs))

    return wrapper


def exists_mention_of_id_type(paper_text: str, regex: str) -> bool:
    pattern = re.compile(regex, re.IGNORECASE)
    id_in_paper = bool(pattern.search(paper_text))
    return id_in_paper


@_standardize
def get_instances_of_id(paper_text: str, regex: str) -> list:
    # Not case sensitive
    pattern = re.compile(regex, flags=re.IGNORECASE | re.VERBOSE)
    instances = pattern.findall(paper_text)
    return instances
