"""
Interface with the protein data bank (PDB) and GenBank.
"""
from functools import cache
from urllib.error import HTTPError
from os import getcwd
import pypdb
from ABDB import database
from Bio.PDB.PDBList import PDBList
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from .anarci_interface import check_if_antibody, extract_VH_VL
from anarci import run_anarci
from anarci import number

PDB_FILEPATH = "./webscraping/proteinids/PDBs/"

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

class GenBankID:
    """Class representating a GenBankID."""

    def __init__(self, genbank_id):
        """Simple initialization, use an Oxford email for
        Entrez."""
        Entrez.email = "vice-chancellor@admin.ox.ac.uk"
        self.genbank_id = genbank_id

    def __bool__(self):
        return bool(self.handle)

    @property
    def citation_info(self):
        """Get the citation info."""
        return self.handle['GBSeq_references'][0]

    @property
    def handle(self):
        """Get the handle."""
        return GenBankID.get_handle(self.genbank_id)

    @property
    def doi(self):
        """Get the doi."""
        doi_info = self.citation_info.get('GBReference_xref', None)
        if not doi_info:
            return
        if doi_info[0]['GBXref_dbname'] == 'doi':
            return doi_info[0]['GBXref_id']

    @property
    def pmid(self):
        """Get the pmid."""
        return self.citation_info.get('GBReference_pubmed', None)

    @property
    def authors(self):
        """Get the authors."""
        return self.citation_info.get('GBReference_authors', None)

    @property
    def sequence(self):
        """Get the sequence. Convert the letters to uppercase."""
        # Note that sequences can be either DNA or amino acid sequences
        # I think dependent on whether it is a nucleotide or protein.
        sequence = self.handle.get('GBSeq_sequence', None)
        if sequence:
            return sequence.upper()

    def get_sequence_strand_amino(self):
        """Get whether it is a heavy or light chain.
        Default scheme is 'imgt'."""
        numbering, chain_type = number(self.sequence, scheme='imgt', allow=set(["H","K","L"]))
        if numbering:
            # replace the Kappa annotation with a light annotation
            # (will come back as L for a lambda chain already).
            chain_type.replace("K","L")
            return chain_type

    @property
    def needs_translation(self):
        """Check if a sequence is in DNA or RNA and thus needs to be
        translated to amino acids.
        We do this by checking if the sequence is a subset of DNA or
        RNA letters."""
        dna_rna_set = {'A', 'C', 'T', 'G', 'U'}
        sequence_set = set(self.sequence)
        return sequence_set.issubset(dna_rna_set)

    def get_protein_ids_non_amino(self):
        """Get the associated protein IDs if these are needed."""
        if not self.needs_translation:  # If already in amino acids, escape
            return
        feature_table = self.handle.get('GBSeq_feature-table', None)
        if not feature_table:
            return
        cds_filter = list(filter(lambda entry: entry['GBFeature_key'] == 'CDS', feature_table))
        protein_ids = list()
        for entry in cds_filter:
            for qual in entry['GBFeature_quals']:
                if qual['GBQualifier_name'] == 'protein_id':
                    protein_ids.append(qual['GBQualifier_value'])
        protein_ids = list(set(protein_ids))
        return [(id_value, 'genbank_protein_id') for id_value in protein_ids]

    @property
    def description(self):
        """Get the description of the GenBank ID"""
        return self.handle.get('GBSeq_definition', None)

    def get_if_antibody_amino(self):
        """Check if the GenBankID is an antibody. Note that this only returns a boolean
        if the sequence is already amino acids.
        We do this by checking if a field in the output from run_anarci is empty where if
        it is empty, it is not an antibody."""
        # if self.needs_translation:  # If not in amino acids, escape
        #     return
        output = run_anarci(self.sequence)
        return bool(output[1][0])

    @staticmethod
    @cache
    def get_handle(genbank_id):
        """
        Get the Entrez information on the genbank_id.
        We allow any entry to be a nucleotide or a protein.
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
