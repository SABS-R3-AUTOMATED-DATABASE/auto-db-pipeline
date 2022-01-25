"""
File containing code for checking against the protein data bank (PDB).
"""
import tempfile
import re
from functools import cache
from Bio.PDB.PDBList import PDBList

TOP_NUM_AUTHORS = 3

@cache
def get_pdb_hash() -> dict:
    """
    Get a hash map of all the existing PDB IDs, such that an existing PDB will 
    return True when looked up in the hash map. This allows for O(1) lookup. 
    We cache this because this can take between 15-20 seconds to load all the existing PDBs.
    """
    pdbl = PDBList(verbose=False, obsolete_pdb="None")  # negligible time
    hash_map = {pdb_id: True for pdb_id in pdbl.get_all_entries()}  # takes a while
    return hash_map


class PDBChecker:
    """
    Object that contains information on the protein data bank (PDB) and interacts with the 
    protein data bank, checking if a PDB ID exists on the data bank, and checking the 
    authors of a particular PDB ID.
    """

    def __init__(self):
        """
        First we store all the existing pdb IDs as a dictionary for O(1) lookup. 
        There are 184,929 IDs as of 2021-12-8 and the retrieval using biopython takes about 7 seconds.
        For some reason, calling PDBList() creates an empty folder in the directory called "obsolete", 
        but this goes away by setting the `obsolte_pdb` parameter to some random string, which I made "None".
        """
        self.pdbl = PDBList(verbose=False, obsolete_pdb="None")
        self.existing_pdbs = get_pdb_hash()

    def get_actual(self, possible_pdbs: list, verbose=True) -> list:
        """
        Takes a list of possible PDB IDs as input. 
        Returns a list of the actual PDB IDs, i.e. the ones from the input list that exist on the PDB database.
        
        Warning: Please remember that html gobble can include actual PDB IDs by chance. So just because a possible
        PDB ID from the paper url html turns out to be an actual PDB ID (is actually on the database), does not 
        mean it was meant to be written in the text of the paper. 
        """
        actual_pdbs = [pdb_id for pdb_id in possible_pdbs if self.existing_pdbs.get(pdb_id, False)]
        if verbose: 
            # Add logger of info or warning of the non-actual PDBs
            print("Out of the", len(possible_pdbs), 'possible PDB IDs scraped', len(actual_pdbs), 'are actual PDB IDs.')
        return actual_pdbs

    @cache
    def get_top_authors(self, pdb_id: str, top_num=TOP_NUM_AUTHORS, verbose=True) -> list:
        """
        Takes an actual PDB ID as input.
        Returns the last names of the top authors for that paper, as retrieved from the PDB database. 
        
        Misc notes:
        If an author has essentially two last names, like "von Kuegelgen", the function will treat those
        as two separate last names. Though this shouldn't matter for practical purposes. 

        The logic of returning only the top authors is that sometimes institutions are named as authors, 
        for example "Seattle Structural Genomics Center for Infectious Disease (SSGCID), McGuire, A.T., Veesler, D."
        or "Midwest Center for Structural Genomics". 
        If there are less than `top_num` authors on the paper, it will return all the authors
        of the paper. 
        """
        temp_dir = tempfile.TemporaryDirectory()
        pdb_file = self.pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=temp_dir.name)
        with open(pdb_file, 'r', encoding='utf-8') as read_file:
            author_txt = ' '.join(filter(lambda line: line.split()[0] == "AUTHOR", read_file.read().splitlines()))
        temp_dir.cleanup()
        top_authors = list(filter(lambda word: len(word) > 1 and word != "AUTHOR", re.findall(r"[\w']+", author_txt)))[:top_num]
        if verbose:
            print("Top Authors scraped from PDB Database:", top_authors)
        return top_authors
    
