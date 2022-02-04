"""
Interface with the protein data bank (PDB).
"""
from functools import cache
import pypdb
from Bio.PDB.PDBList import PDBList

class PdbID:
    """
    Class of type pdb.
    """
    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        self.citation_info = PdbID.get_citation_info(self.pdb_id)

    @staticmethod
    @cache
    def get_citation_info(pdb_id: str) -> dict:
        return pypdb.get_info(pdb_id)['citation'][0]

    @property
    def doi(self) -> str:
        return self.citation_info['pdbx_database_id_doi']

    @property
    def pmid(self) -> str:
        return self.citation_info['pdbx_database_id_pub_med']

    @property
    def authors(self) -> list:
        return self.citation_info['rcsb_authors']

    @property
    def exists(self) -> bool:
        """
        Check whether the possible `pdb_id` exists on the database,
        (i.e. is an actual `pdb_id`).
        """
        existing_pdbs = PdbID.get_pdb_hash()
        return existing_pdbs.get(self.pdb_id, False)


    @staticmethod
    @cache
    def get_pdb_hash() -> dict:
        """
        Get a hash map of all the existing PDB IDs, such that an existing PDB will
        return True when looked up in the hash map. This allows for O(1) lookup.
        We cache this because this can take between 15-20 seconds to
        load all the existing PDBs.
        """
        pdbl = PdbID._get_PDBList()
        # Takes a while to load
        hash_map = {pdb_id: True for pdb_id in pdbl.get_all_entries()}
        return hash_map

    @staticmethod
    def _get_PDBList():
        """
        Gets a PDBList object (in a controlled way).
        For some reason, calling PDBList() creates an empty folder
        in the directory called "obsolete",
        but this goes away by setting the `obsolte_pdb` parameter to
        some random string, which I made "None".

        This takes negligible time to generate.
        """
        return PDBList(verbose=False, obsolete_pdb="None")

