"""
Interface with the protein data bank (PDB).
"""
from functools import cache
import pypdb
from ABDB import database
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from .anarci_interface import check_if_antibody, extract_VH_VL

PDB_FILEPATH = "./protein_scrape/PDBs/"

class PdbID:
    """Class representating a PdbID."""

    def __init__(self, pdb_id):
        """Simple initialization, store the pdb_id value."""
        self.pdb_id = pdb_id

    def __bool__(self) -> bool:
        """Check whether the possible `pdb_id` exists on the database."""
        existing_pdbs = PdbID.get_pdbs()
        return self.pdb_id in existing_pdbs

    @property
    def pdb_file(self):
        """Get the pdb file.
        This automatically caches if the file already exists in the filepath."""
        pdb_list = PdbID.get_pdb_list()
        return pdb_list.retrieve_pdb_file(self.pdb_id, file_format="pdb", pdir=PDB_FILEPATH)

    @property
    def sequence_repr(self):
        """Get the sequence representation of a PDB."""
        parse = SeqIO.parse(self.pdb_file, 'pdb-atom')
        seq_repr = {str(record.id): str(record.seq).replace('X', '') for record in parse}
        return seq_repr

    @property
    def is_antibody(self) -> bool:
        """Implement the multiple ways to check if the PDB is an antibody."""
        if self.on_sabdab:
            return True
        if self.check_anarci_if_antibody():
            return True
        return False

    @property
    def vh_vl_anarci(self) -> dict:
        """Get the VH and VL sequences. This will only work if it's an antibody."""
        return extract_VH_VL(self.sequence_repr)

    @cache
    def check_anarci_if_antibody(self) -> bool:
        """Check with ANARCI if the PDB is an antibody. Note the asssumption
        sated in the function's docstring."""
        return check_if_antibody(self.sequence_repr)

    @property
    def handle(self):
        """Get the handle."""
        return PdbID.get_handle(self.pdb_id)

    @property
    def citation_info(self) -> dict:
        """Get the citation info."""
        return self.handle['citation'][0]

    @property
    def doi(self) -> str:
        """Get the doi."""
        return self.citation_info.get('pdbx_database_id_doi', None)

    @property
    def pmid(self) -> str:
        """Get the pmid."""
        return self.citation_info.get('pdbx_database_id_pub_med', None)

    @property
    def authors(self) -> list:
        """Get the authors."""
        return self.citation_info['rcsb_authors']

    @staticmethod
    @cache
    def get_handle(pdb_id):
        """Retrieve the handle and cache it (takes time)."""
        return pypdb.get_info(pdb_id)

    @property
    def on_sabdab(self) -> bool:
        """Returns whether the pdb_id is on SABDAB."""
        return bool(database.fetch(self.pdb_id))

    @staticmethod
    @cache
    def get_pdbs() -> set:
        """
        Get a set of all PDB IDs. This allows for O(1) lookup.
        We cache this because this can take between 15-20 seconds to
        load all the existing PDBs.

        For some reason, calling PDBList() creates an empty folder
        in the directory called "obsolete",
        but this goes away by setting the `obsolte_pdb` parameter to
        some random string, which I made "None".

        """
        pdbl = PdbID.get_pdb_list()
        # Takes a while to load
        all_entries = pdbl.get_all_entries()  # list
        return set(all_entries)

    @staticmethod
    def get_pdb_list():
        """Get a PDBList object in a controlled way."""
        return PDBList(verbose=False, obsolete_pdb="None")
