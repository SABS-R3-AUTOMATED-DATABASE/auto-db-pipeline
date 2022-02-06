"""
Interface with the protein data bank (PDB) and GenBank.
"""
from functools import cache
import pypdb
from urllib.error import HTTPError
from Bio.PDB.PDBList import PDBList
from Bio import Entrez

class PdbID:
    """Class representating a PdbID."""

    def __init__(self, pdb_id):
        self.pdb_id = pdb_id

    def __bool__(self) -> bool:
        """Check whether the possible `pdb_id` exists on the database."""
        existing_pdbs = PdbID.get_pdbs()
        return self.pdb_id in existing_pdbs

    @property
    def citation_info(self):
        return PdbID.get_citation_info(self.pdb_id)

    @staticmethod
    @cache
    def get_citation_info(pdb_id: str) -> dict:
        return pypdb.get_info(pdb_id)['citation'][0]

    @property
    def doi(self) -> str:
        return self.citation_info.get('pdbx_database_id_doi', None)

    @property
    def pmid(self) -> str:
        return self.citation_info.get('pdbx_database_id_pub_med', None)

    @property
    def authors(self) -> list:
        return self.citation_info['rcsb_authors']

    @property
    def sequence(self):
        # add seqeuence functionality
        pass


    @staticmethod
    @cache
    def get_pdbs() -> set:
        """
        Get a set of all PDB IDs. This allows for O(1) lookup.
        We cache this because this can take between 15-20 seconds to
        load all the existing PDBs.
        """
        pdbl = PdbID._get_PDBList()
        # Takes a while to load
        return set(pdbl.get_all_entries())

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


class GenBankID:
    """Class representating a GenBankID."""

    def __init__(self, genbank_id):
        Entrez.email = "vice-chancellor@admin.ox.ac.uk"
        self.genbank_id = genbank_id

    def __bool__(self):
        return bool(self.handle)

    @property
    def citation_info(self):
        return self.handle['GBSeq_references'][0]

    @property
    def handle(self):
        return GenBankID.get_handle(self.genbank_id)

    @property
    def doi(self):
        doi_info = self.citation_info.get('GBReference_xref', None)
        if not doi_info:
            return None
        if doi_info[0]['GBXref_dbname'] == 'doi':
            return doi_info[0]['GBXref_id']

    @property
    def pmid(self):
        return self.citation_info.get('GBReference_pubmed', None)

    @property
    def authors(self):
        return self.citation_info.get('GBReference_authors', None)

    @property
    def sequence(self):
        return self.handle.get('GBSeq_sequence', None)

    @staticmethod
    @cache
    def get_handle(genbank_id):
        """
        Get the Entrez information on the genbank_id.
        """
        try:
            handle = Entrez.efetch(id=genbank_id, db='protein',
                                    rettype="gb", retmode="xml")
        except HTTPError:
            try:
                handle = Entrez.efetch(id=genbank_id, db='nucleotide',
                                     rettype="gb", retmode="xml")
            except HTTPError:
                # e.g.  incorrect RefSeq type or does not exist
                return None

        return Entrez.read(handle)[0]
